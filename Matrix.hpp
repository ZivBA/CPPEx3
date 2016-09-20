//
// Created by rooty on 9/19/16.
//

#ifndef CPPEX3_Matrix_HPP
#define CPPEX3_Matrix_HPP

#include <vector>
#include <stdexcept>
#include <thread>
#include <future>
#include "Complex.h"

// undefine to disable range checking
#define RANGE_CHECK

template<class T>
class Matrix
{

public:
    unsigned myRows;  // number of myRows
    unsigned myCols;  // number of columns

private:

    static bool isParallel;
    std::vector<T> elements; // array of elements
    // number of available logical threads;
    unsigned concurentThreadsSupported = std::thread::hardware_concurrency();


protected:
    // range check function for Matrix access
    void range_check(unsigned i, unsigned j) const;

public:

    static void setParallel(bool isPar)
    {
        if (isPar != isParallel)
        {
            Matrix::isParallel = isPar;
            std::cout <<
                      "Generic Matrix mode changed to " <<
                      (isPar ? "Parallel" : "non-Parallel") <<
                      " mode." << std::endl;
        }
    }

    typedef typename std::vector<T>::const_iterator const_iterator;


    T &operator()(unsigned i, unsigned j)
    {
#ifdef RANGE_CHECK
        range_check(i, j);
#endif
        return elements[i * myCols + j];
    }

    const T &operator()(unsigned i, unsigned j) const
    {
#ifdef RANGE_CHECK
        range_check(i, j);
#endif
        return elements[i * myCols + j];
    }

    const T &element(unsigned i, unsigned j) const
    {
#ifdef RANGE_CHECK
        range_check(i, j);
#endif
        return elements[i * myCols + j];
    }

    T &element(unsigned i, unsigned j)
    {
#ifdef RANGE_CHECK
        range_check(i, j);
#endif
        return elements[i * myCols + j];
    }

    static std::vector<T>
    multHelper(const Matrix<T> &pMatrix, const Matrix<T> &pMatrix1, const int i, const int m,
               const int p) ;

    static std::vector<T>
    addHelper(const Matrix<T> &a, const Matrix<T> &b, const int i, const int m) ;

    // constructors
    Matrix();

    Matrix(unsigned rows, unsigned cols);

    Matrix(const Matrix<T> &);

    Matrix( Matrix<T> &&);

    Matrix(unsigned int rows, unsigned int cols, const std::vector<T> &cells);

    // destructor
    ~Matrix();

    // assignment
    Matrix<T> &operator=(const Matrix<T> &);

    // comparison
    bool operator==(const Matrix<T> &) const;

    bool operator!=(const Matrix<T> &oMat) const
    {
        return !(*this == oMat);
    }

    // addition/subtraction

    Matrix<T> operator+(const Matrix<T> &oMat) const;

    Matrix<T> operator-(const Matrix<T> &oMat) const;

    // Matrix multiplication
    Matrix<T> operator*(const Matrix<T> &) const;

    // Other stuff
    bool isSquareMatrix() const;

     unsigned int rows() const
    {
        return myRows;
    }

     unsigned int cols()
    {
        return myCols;
    }


    Matrix<T> trans() const;
    //iterator
    const_iterator begin();

    const_iterator end();

    //stream
    friend std::ostream& operator<<(std::ostream& os, const Matrix<T>& mat)
    {
        for (unsigned i = 0; i < mat.myRows; i++)
        {
            for (unsigned j = 0; j < mat.myCols; j++)
            {
                os << mat(i,j) << "\t";
            }
            os << std::endl;
        }
        return os;
    }
};


template<class T>
Matrix<T>::Matrix() :  myRows(1) ,myCols(1), elements((0.0))
{
}

template<class T>
Matrix<T>::Matrix(unsigned rows, unsigned cols)
        : myRows(rows), myCols(cols), elements(rows * cols, T(0.0))
{
    if (rows == 0 || cols == 0)
    {
        throw std::range_error("attempt to create a degenerate Matrix");
    }
};

template<class T>
Matrix<T>::Matrix(const Matrix<T> &cp)
        : myRows(cp.myRows), myCols(cp.myCols), elements(cp.elements)
{
}

template<class T>
Matrix<T>::Matrix(unsigned int rows, unsigned int cols, const std::vector<T> &cells) :
        myRows(rows), myCols(cols), elements(cells)
{

}

template<class T>
Matrix<T>::Matrix(Matrix &&mv) : myRows(mv.myRows), myCols(mv.myCols)
{
    using std::swap;
    std::swap(this->elements, mv.elements);

}

template<class T>
Matrix<T>::~Matrix()
{
}

template<class T>
Matrix<T> &Matrix<T>::operator=(const Matrix<T> &cp)
{
    this->myCols = cp.myCols;
    this->myRows = cp.myRows;
    elements.clear();
    elements.assign(cp.elements.begin(),cp.elements.end());

    return *this;
}

template<class T>
void Matrix<T>::range_check(unsigned i, unsigned j) const
{
    if (myRows <= i)
    {
        throw std::range_error("Matrix access row out of range");
    }
    if (myCols <= j)
    {
        throw std::range_error("Matrix access col out of range");
    }
}
template<class T>
bool Matrix<T>::isSquareMatrix() const
{
    return myCols == myRows;
}

template<class T>
bool Matrix<T>::operator==(const Matrix &oMat) const
{
    return this->elements == oMat.elements;


}

template<class T>
Matrix<T> Matrix<T>::operator*(const Matrix<T> &oMat) const
{
    if (this->myCols != oMat.myRows)
    {
        throw std::range_error("Matrix dimensions incompatible");
    }

    int n = this->myRows;
    int m = this->myCols;
    int p = oMat.myCols;
    Matrix<T> temp = Matrix(n, p);

    int numthreads = isParallel ? concurentThreadsSupported : 1;
    int rowsToProc = n;
    std::vector<std::future<std::vector<T>>> threads;

    temp.elements.clear();

    int i=0;
    while (rowsToProc > 0)
    {
        for (int j =i; j< i+numthreads && j < n; j++)
        {
            threads.push_back(std::async(Matrix<T>::multHelper,*this, oMat, j, m, p));
        }

        for (int j =i; j< i+numthreads && j < n; j++)
        {
            std::vector<T> tempVec = threads[j].get();
            temp.elements.insert(temp.elements.end(), tempVec.begin(), tempVec.end());
        }
        i+= numthreads;
        rowsToProc -= numthreads;
    }

//    for (int i = 0; i < n; i++)
//    {
//        threads.push_back(std::async(Matrix<T>::multHelper,*this, oMat, i, m, p));
//    }
//    temp.elements.clear();
//    for (int i = 0; i < n; i++)
//    {
//        std::vector<T> tempVec = threads[i].get();
//        temp.elements.insert(temp.elements.end(), tempVec.begin(), tempVec.end());
//    }

    return temp;


}


template<class T>
 std::vector<T>
Matrix<T>::multHelper(const Matrix &a, const Matrix &b, const int i, const int m, const int p)
{
    std::vector<T> resultLine;
    for (int j = 0; j < p; j++)
    {
        T sum = 0.0;
        for (int k = 0; k < m; k++)
        {
            sum += a(i, k) * b(k, j);
        }
        resultLine.push_back(sum);
    }
    return resultLine;

}

template<class T>
Matrix<T> Matrix<T>::operator+(const Matrix<T> &oMat) const
{
    if (this->myCols != oMat.myCols || this->myRows != oMat.myRows)
    {
        throw std::range_error("Cannot add matrices of different dimensions.");
    }
    int n = this->myRows;
    int m = this->myCols;
    int numthreads = isParallel ? concurentThreadsSupported : 1;
    int rowsToProc = n;
    std::vector<std::future<std::vector<T>>> threads;

    Matrix<T> temp = Matrix(n, m);
    temp.elements.clear();

    int i=0;
    while (rowsToProc > 0)
    {
        for (int j =i; j< i+numthreads && j < n; j++)
        {
            threads.push_back(std::async(Matrix<T>::addHelper,*this, oMat, j, m));
        }

        for (int j =i; j< i+numthreads && j < n; j++)
        {
            std::vector<T> tempVec = threads[j].get();
            temp.elements.insert(temp.elements.end(), tempVec.begin(), tempVec.end());
        }
        i+= numthreads;
        rowsToProc -= numthreads;
    }


//    for (int i = 0; i < n; i++)
//    {
//        threads.push_back(std::async(Matrix<T>::addHelper,*this, oMat, i, m));
//    }
//
//    for (int i = 0; i < n; i++)
//    {
//        std::vector<T> tempVec = threads[i].get();
//        temp.elements.insert(temp.elements.end(), tempVec.begin(), tempVec.end());
//    }

    return temp;


}

template<class T>
 std::vector<T>
Matrix<T>::addHelper(const Matrix &a, const Matrix &b, const int i, const int m)

{
    std::vector<T> tempVec;
    for (int j = 0; j < m; j++)
    {
        tempVec.push_back(a(i, j) + b(i, j));
    }
    return tempVec;
}

template<class T>
Matrix<T> Matrix<T>::operator-(const Matrix &oMat) const
{
    if (this->myCols != oMat.myCols || this->myRows != oMat.myRows)
    {
        throw std::range_error("Cannot add matrices of different dimensions.");
    }
    Matrix<T> temp(this->myRows, this->myCols);

    for (int i = 0; i < this->myRows; i++)
    {
        for (int j = 0; j < this->myCols; j++)
        {
            temp.elements[i * myCols + j] = this->operator()(i, j) - oMat(i, j);
        }
    }
    return temp;
}


template<class T>
Matrix<T> Matrix<T>::trans() const
{
    int n = this->myRows;
    int m = this->myCols;

    std::vector<T> newMatVec;
    for (int i = this->myRows - 1; i >= 0; i--)
    {
        for (int j = this->myCols-1; j >= 0; j--)
        {
            newMatVec.push_back(elements[i * myCols + j]);
        }
    }

    return Matrix(n, m, newMatVec);

}

template<class T>
typename std::vector<T>::const_iterator Matrix<T>::begin()
{
    return elements.begin();
}

template<class T>
typename std::vector<T>::const_iterator Matrix<T>::end()
{
    return elements.end();
}



template<> bool Matrix<Complex>::isParallel = false;

template<>
Matrix<Complex> Matrix<Complex>::trans() const
{

    int n = this->myRows;
    int m = this->myCols;

    std::vector<Complex> newMatVec;
    for (int i = this->myRows - 1; i >= 0; i--)
    {
        for (int j = this->myCols-1; j >= 0; j--)
        {
            newMatVec.push_back(elements[i * myCols + j]);
        }
    }



    Matrix<Complex> tempMat(n, m, newMatVec);


    for (unsigned i = 0; i < this->myRows; i++)
    {
        for (unsigned j = 0; j < this->myCols; j++)
        {
            tempMat(i,j) = tempMat(i,j).conj() ;
        }
    }
    return tempMat;

}

#endif //CPPEX3_Matrix_HPP
