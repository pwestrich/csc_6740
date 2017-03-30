
#include <cmath>
#include <cstdlib>
#include <cstdint>
#include <fstream>
#include <iostream>
#include <string>
#include <utility>

#include <mpi.h>

#include "matrix.h"

#define MASTER 0 //the master process is 0

/*
*	Holds the parameters of the problem so my argument lists aren't long.
*/
struct Parameters {

	uint32_t rank;
	uint32_t p;
	uint32_t n;
	uint32_t blockSize;

};

/*
*	Prints the usage of the program given its name from argv[0].
*/
inline void usage(const char *name){

	std::cout << "Usage: mpiexec -n <num_procs> " << name << " <matrix_size>" << std::endl;
	std::cout << "where <matrix_size> % sqrt(<num_procs>) = 0)" << std::endl;
	std::cout << "  and <matrix_size> >= 1" << std::endl;

}

/*
*	Multiplies two matricies such that A * B = C, using Cannon's algorithm.
*/
void mpi_cannon_multiply(const Parameters &params){

	//Step 1: Partition the data
	Matrix<double> *myA = Matrix<double>::initializeMatrix(params.blockSize);
	Matrix<double> *myB = Matrix<double>::initializeMatrix(params.blockSize);
	Matrix<double> *myC = new Matrix<double>(params.blockSize);

	Matrix<double> *nextA = new Matrix<double>(params.blockSize);
	Matrix<double> *nextB = new Matrix<double>(params.blockSize);

	//Step 2: Split processes into their rows and columns
	const uint32_t rowSize = static_cast<uint32_t>(sqrt(params.p));
	int rowRank, colRank;

	MPI_Comm rowComm, colComm;
	MPI_Comm_split(MPI_COMM_WORLD, params.rank / rowSize, params.rank, &rowComm);
	MPI_Comm_split(MPI_COMM_WORLD, params.rank % rowSize, params.rank, &colComm);
	MPI_Comm_rank(rowComm, &rowRank);
	MPI_Comm_rank(colComm, &colRank);

	std::ofstream outFile(std::to_string(params.rank) + ".txt");
	outFile << "rank: " << params.rank << "("<< rowRank << "," << colRank << ")"<< std::endl;
	outFile << "A ";
	myA->print(outFile);
	outFile << "B ";
	myB->print(outFile);

	const uint32_t nextRow = (rowRank + 1) % rowSize;
	const uint32_t nextCol = (colRank + 1) % rowSize;
	const uint32_t prevRow = (rowRank - 1) % rowSize;
	const uint32_t prevCol = (colRank - 1) % rowSize;

	//we need to do these once per block
	for (uint32_t block = 0; block < rowSize; ++block){

		//Step 3: Start getting the next set
		MPI_Request requests[4];
		MPI_Status statuses[4];
		MPI_Isend(static_cast<void*>(myA->getRaw()), myA->getSize(), MPI_DOUBLE, nextRow, 0, rowComm, &requests[0]);
		MPI_Isend(static_cast<void*>(myB->getRaw()), myB->getSize(), MPI_DOUBLE, nextCol, 0, colComm, &requests[1]);

		MPI_Irecv(static_cast<void*>(nextA->getRaw()), nextA->getSize(), MPI_DOUBLE, prevRow, 0, rowComm, &requests[2]);
		MPI_Irecv(static_cast<void*>(nextB->getRaw()), nextB->getSize(), MPI_DOUBLE, prevCol, 0, colComm, &requests[3]);

		//Step 4: Multiply my set together (naive method for now)
		for (uint32_t i = 0; i < params.blockSize; ++i){

			for (uint32_t j = 0; j < params.blockSize; ++j){

				double sum = 0.0;
				for (uint32_t k = 0; k < params.blockSize; ++k){

					sum += myA->getElement(i, k) * myB->getElement(k, j);

				}

				myC->addToElement(i, j, sum);

			}

		}

		//Step 5: Wait on step 3 to finish, then swap and loop
		MPI_Waitall(4, requests, statuses);

		std::swap(myA, nextA);
		std::swap(myB, nextB);

	}

	//Step 6: write results to file
	outFile << "C ";
	myC->print(outFile);
	outFile.close();

	//Step 7: Clean up.
	MPI_Comm_free(&rowComm);
	MPI_Comm_free(&colComm);
	delete myA;
	delete myB;
	delete myC;
	delete nextA;
	delete nextB;

}

int main(int argc, char *argv[]){

	//initialize MPI, find out my rank and size
    MPI_Init(&argc, &argv);
	int _rank = 0;
	int _p 	= 0;
	MPI_Comm_size(MPI_COMM_WORLD, &_p);
	MPI_Comm_rank(MPI_COMM_WORLD, &_rank);

	if (argc < 2){

		if (_rank == MASTER){

			std::cerr << "Error: Too few arguments." << std::endl;
			usage(argv[0]);

		}

		MPI_Finalize();
		return EXIT_FAILURE;

	}

	const uint32_t _n = std::stoul(argv[1]);
	const Parameters params = {static_cast<uint32_t>(_rank), static_cast<uint32_t>(_p), _n, _n / _p};

	//test arguments
	if (params.n <= 0){

		if (params.rank == MASTER){

			std::cerr << "Error: n is negative." << std::endl;
			usage(argv[0]);

		}

		MPI_Finalize();
		return EXIT_FAILURE;

	} else if (sqrt(params.p) != floor(sqrt(params.p))){

		if (params.rank == MASTER){

			std::cerr << "Error: sqrt(" << params.p << ") is not an even number." << std::endl;
			usage(argv[0]);

		}

		MPI_Finalize();
		return EXIT_FAILURE;

	} else if ((params.n % static_cast<uint32_t>(floor(sqrt(params.p)))) != 0){

		if (params.rank == MASTER){

			std::cerr << "Error: " << params.n << " % sqrt(" << params.p << ") != 0" << std::endl;
			usage(argv[0]);

		}

		MPI_Finalize();
		return EXIT_FAILURE;

	} else if (params.rank == MASTER) {

		std::cout << "Multiplying two dense matricies of size " << params.n << " with block size "
				  << params.blockSize << " and " << params.p << " processes." << std::endl;

	}

	mpi_cannon_multiply(params);

	//exit
	MPI_Finalize();
	return EXIT_SUCCESS;

}
