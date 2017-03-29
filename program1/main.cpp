
#include <cmath>
#include <cstdlib>
#include <cstdint>
#include <iostream>

#include <mpi.h>

#include "matrix.h"

#define MASTER 0 //the master process is 0

inline void usage(const char *name){

	std::cout << "Usage: mpiexec -n <num_procs> " << name << " <matrix_size>" << std::endl;
	std::cout << "where <matrix_size> % sqrt(<num_procs>) = 0)" << std::endl;
	std::cout << "  and <matrix_size> >= 1" << std::endl;

}

void mpi_multiply(const uint32_t rank, const uint32_t n, const uint32_t blockSize, const uint32_t p){

	//initialize matricies
	Matrix<double> *myA = Matrix<double>::identityMatrix(blockSize);
	Matrix<double> *myB = Matrix<double>::identityMatrix(blockSize);
	Matrix<double> *myC = new Matrix<double>(blockSize);

	//test the single threaded multiply
	Matrix<double>::multiply(*myA, *myB, *myC);
	assert(*myC == *myA);
	assert(*myC == *myB);

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

	const uint32_t rank 		= _rank;
	const uint32_t p 			= _p;
	const uint32_t n 			= std::stoul(argv[1]);
	const uint32_t blockSize  	= n / p;

	//test arguments
	if (n <= 0){

		if (rank == MASTER){

			std::cerr << "Error: n is negative." << std::endl;
			usage(argv[0]);

		}

		MPI_Finalize();
		return EXIT_FAILURE;

	} else if (sqrt(p) != floor(sqrt(p))){

		if (rank == MASTER){

			std::cerr << "Error: sqrt(" << p << ") is not an even number." << std::endl;
			usage(argv[0]);

		}

		MPI_Finalize();
		return EXIT_FAILURE;

	} else if ((n % static_cast<uint32_t>(floor(sqrt(p)))) != 0){

		if (rank == MASTER){

			std::cerr << "Error: " << n << " % sqrt(" << p << ") != 0" << std::endl;
			usage(argv[0]);

		}

		MPI_Finalize();
		return EXIT_FAILURE;

	} else if (rank == MASTER) {

		std::cout << "Multiplying two dense matricies of size " << n << " with block size "
				  << blockSize << " and " << p << " processes." << std::endl;

	}

	mpi_multiply(rank, n, blockSize, p);

	//exit
	MPI_Finalize();
	return EXIT_SUCCESS;

}
