
#ifndef MATRIX__H
#define MATRIX__H

#include <cassert>
#include <cstdint>
#include <iterator>
#include <random>
#include <vector>

template <class T>
class Matrix {

private:

    uint32_t n;         //#rows/cols
    uint32_t size;      //how many elements?
    std::vector<T> elements; //container holding them

public:

    //Constructs an empty matrix with enough size for n*n elements
    Matrix(const uint32_t _n) : n(_n), size(n * n), elements(size, static_cast<T>(0)){};

    //Fills the matrix with the specified elements
    Matrix(std::input_iterator_tag first, const uint32_t _n) : n(_n), size(n * n), elements(first, first + size){};

    bool operator==(const Matrix<T> &other){

        if (n != other.getN()) return false;

        for (uint32_t row = 0; row < n; ++row){

            for (uint32_t col = 0; col < n; ++col){

                if (getElement(row, col) != other.getElement(row, col)) return false;

            }

        }

        return true;

    }

    static Matrix<T> *identityMatrix(const uint32_t blockSize){

        Matrix<T> *matrix = new Matrix<double>(blockSize);

        for (uint32_t i = 0; i < blockSize; ++i){

            matrix->setElement(i, i, static_cast<T>(1));

        }

        return matrix;

    }

    static Matrix<T> *initializeMatrix(const uint32_t blockSize){

    	std::default_random_engine generator;
    	std::uniform_int_distribution<T> distribution(-10, 10);

    	Matrix<double> *matrix = new Matrix<double>(blockSize);

    	for (uint32_t row = 0; row < blockSize; ++row){

    		for (uint32_t col = 0; col < blockSize; ++col){

    			matrix->setElement(row, col, distribution(generator));

    		}

    	}

    	return matrix;

    }

    //Multples C = A * B using the naive method, single-core
    static void multiply(const Matrix<T> &A, const Matrix<T> &B, Matrix<T> &C){

        const uint32_t n = A.getN();

        assert(n == B.getN());
        assert(n == C.getN());

        for (uint32_t i = 0; i < n; ++i){

            for (uint32_t j = 0; j < n; ++j){

                T sum = static_cast<T>(0);
                for (uint32_t k = 0; k < n; ++k){

                    sum += A.getElement(i, k) * B.getElement(k, j);

                }

                C.setElement(i, j, sum);

            }

        }

    }

    uint32_t getN() const {

        return n;

    }

    uint32_t getSize() const {

        return n * n;

    }

    T getElement(const uint32_t row, const uint32_t col) const {

        return elements[(row * n) + col];

    }

    void setElement(const uint32_t row, const uint32_t col, const T element){

        elements[(row * n) + col] = element;

    }

    T* getRaw() const {

        return elements.data();

    }

};

#endif
