
#ifndef MATRIX__H
#define MATRIX__H

#include <algorithm>
#include <cassert>
#include <cstdint>
#include <iostream>
#include <iterator>
#include <random>
#include <vector>

template <class T>
class Matrix {

private:

    const uint32_t n;       //#rows/cols
    const uint32_t size;    //how many elements?
    T *const elements;      //container holding them

public:

    //Constructs an empty matrix with enough size for n*n elements
    Matrix(const uint32_t _n) : n(_n), size(n * n), elements(new T[size]){

        std::fill(elements, elements + size, 0);

    };

    Matrix(const Matrix& src) = delete;
    Matrix& operator=(const Matrix&) = delete;

    ~Matrix(){

        delete [] elements;

    }

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

    	static std::default_random_engine generator;
    	static std::uniform_int_distribution<int> distribution(-10, 10);

    	Matrix<double> *matrix = new Matrix<double>(blockSize);

    	for (uint32_t row = 0; row < blockSize; ++row){

    		for (uint32_t col = 0; col < blockSize; ++col){

    			matrix->setElement(row, col, distribution(generator));

    		}

    	}

    	return matrix;

    }

    uint32_t getN() const {

        return n;

    }

    uint32_t getSize() const {

        return size;

    }

    T getElement(const uint32_t row, const uint32_t col) const {

        return elements[(row * n) + col];

    }

    void setElement(const uint32_t row, const uint32_t col, const T element){

        assert(row < n);
        assert(col < n);
        elements[(row * n) + col] = element;

    }

    void addToElement(const uint32_t row, const uint32_t col, const T element){

        assert(row < n);
        assert(col < n);
        elements[(row * n) + col] += element;

    }

    T* getRaw(){

        return elements;

    }

    void print(std::ostream &out) const {

        out << n << std::endl;

        for (uint32_t i = 0; i < n; ++i){

            for (uint32_t j = 0; j < n; ++j){

                out << getElement(i, j) << ",";

            }

            out << std::endl;

        }

        out << std::endl;

    }

};

#endif
