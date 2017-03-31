
#include <algorithm>
#include <cassert>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <random>
#include <vector>

#define BSP_DISABLE_LEGACY
#include "bsp/bsp.h"

typedef std::vector<double> Matrix;

inline double getElement(const Matrix &matrix, const uint32_t row, const uint32_t col, const uint32_t n){

    return matrix[(row * n) + col];

}

inline void setElement(Matrix &matrix, const uint32_t row, const uint32_t col, const uint32_t n, const double val){

    matrix[(row * n) + col] = val;

}

inline void addToElement(Matrix &matrix, const uint32_t row, const uint32_t col, const uint32_t n, const double val){

    matrix[(row * n) + col] += val;

}

Matrix randomMatrix(const uint32_t size){

    static std::default_random_engine generator;
    static std::uniform_int_distribution<int> distribution(-10, 10);
    Matrix matrix(size);

    for (uint32_t i = 0; i < size; ++i){

        matrix[i] = distribution(generator);

    }

    return matrix;

}

void writeToFile(std::ostream &out, const Matrix &matrix, const uint32_t n){

    out << n << std::endl;

    for (uint32_t i = 0; i < n; ++i){

        for (uint32_t j = 0; j < n; ++j){

            out << getElement(matrix, i, j, n) << ",";

        }

        out << std::endl;

    }

    out << std::endl;

}

/*
*	Prints the usage of the program given its name from argv[0].
*/
inline void usage(const char *name){

	std::cout << "Usage: ./" << name << " <num_procs> <matrix_size>" << std::endl;
	std::cout << "where <matrix_size> % sqrt(<num_procs>) = 0)" << std::endl;
	std::cout << "  and <matrix_size> >= 1" << std::endl;

}

int main(const int argc, const char *argv[]){

    //test arguments
    if (argc < 3){

        std::cerr << "Too few arguments." << std::endl;
        usage(argv[0]);
        return EXIT_FAILURE;

    }

    const uint32_t p                = std::stoul(argv[1]);
    const uint32_t n                = std::stoul(argv[2]);
    const uint32_t blocksPerRow     = static_cast<uint32_t>(sqrt(p));
    const uint32_t blockSize        = n / p;
    const uint32_t blockMatrixSize  = blockSize * blockSize;

    if (n <= 0){

        std::cerr << "Error: n is negative." << std::endl;
        usage(argv[0]);

    } else if (sqrt(p) != floor(sqrt(p))){

        std::cerr << "Error: sqrt(" << p << ") is not an even number." << std::endl;
        usage(argv[0]);

    } else if ((n % static_cast<uint32_t>(floor(sqrt(p)))) != 0){

        std::cerr << "Error: " << n << " % sqrt(" << p << ") != 0" << std::endl;
        usage(argv[0]);

    } else {

        std::cout << "Multiplying two dense matricies of size " << n << " with block size "
                  << blockSize << " and " << p << " processes." << std::endl;

    }

    BSPLib::Execute([&]{

        //Step 1: Calculate my indicies
        const uint32_t pid      = BSPLib::ProcId();
        const uint32_t row      = pid / blocksPerRow;
        const uint32_t leftPid  = ((pid - 1) % blocksPerRow) + (row * blocksPerRow);
        const uint32_t rightPid = ((pid + 1) % blocksPerRow) + (row * blocksPerRow);
        const uint32_t upPid    = (pid - blocksPerRow) % p;
        const uint32_t downPid  = (pid + blocksPerRow) % p;

        //Step 2: Reserve some space
        Matrix myA = randomMatrix(blockMatrixSize);
        Matrix myB = randomMatrix(blockMatrixSize);
        Matrix myC(blockMatrixSize, 0);
        Matrix nextA(blockMatrixSize);
        Matrix nextB(blockMatrixSize);

        std::ofstream outFile(std::to_string(pid) + ".txt");
        outFile << "rank: " << pid << std::endl;
        outFile << "A ";
        writeToFile(outFile, myA, blockSize);
        outFile << "B ";
        writeToFile(outFile, myB, blockSize);

        BSPLib::PushContainer(myA);
        BSPLib::PushContainer(myB);
        BSPLib::PushContainer(nextA);
        BSPLib::PushContainer(nextB);
        BSPLib::Sync();

        for (uint32_t block = 0; block < blocksPerRow; ++block){

            //Step 3: Multiply my matricies
            for (uint32_t i = 0; i < blockSize; ++i){

                for (uint32_t j = 0; j < blockSize; ++j){

                    double sum = 0.0;
                    for (uint32_t k = 0; k < blockSize; ++k){

                        sum += getElement(myA, i, k, blockSize) * getElement(myB, k, j, blockSize);

                    }

                    addToElement(myC, i, j, blockSize, sum);

                }

            }

            BSPLib::Sync();

            //Step 4: Shift
            BSPLib::PutContainer(rightPid, myA);
            BSPLib::PutContainer(downPid, myB);
            BSPLib::GetContainer(leftPid, nextA);
            BSPLib::GetContainer(upPid, nextB);
            BSPLib::Sync();

            //Step 5: Swap
            std::swap(myA, nextA);
            std::swap(myB, nextB);

        }

        //Step 5: Clean up and exit
        outFile << "C ";
        writeToFile(outFile, myC, blockSize);
        outFile.close();


    }, p);
    return 0;

}
