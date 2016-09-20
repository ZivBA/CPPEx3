/**
 * A driver for testing the implementation of the Matrix.hpp file which is
 * an implementationf of a generic (template) Matrix
 */
#include <iostream>
#include <cassert>
#include <vector>
#include <string>
#include <cstdlib>
#include <sstream>
#include <iostream>
#include <sstream>


#include "Matrix.hpp"
#include "Complex.h"

#define LINE "=========="
#define MAT_LINE "--------"
#define DELIM ','

#define ADD 1
#define MUL 2
#define TRANSPOSE 3

#define INT 1
#define DOUBLE 2
#define COMPLEX 3

#define NUM_KNOWN_FIELDS 3

const std::string OPERATIONS[6] = {"Matrix + Matrix", "Matrix * Matrix", "transpose"};
const int NUM_OPS = 3;

const std::string FIELD_NAMES[3] = {"ints", "Doubles", "Complexs"};

std::string g_line;

template <typename T>
void unaryOperation(const int operatorChoice);

template <typename T>
void binaryOperation(const int operatorChoice);

template <typename T>
void printResultMatrix(const Matrix<T>& mat);

template <typename T>
void readMatrixInfo(int& rows, int& cols, std::vector<T>& cells);

template <typename T>
void getNumFromString(const std::string &str, T *num);

template <>
void getNumFromString(const std::string &str, int *num);

template <typename T>
const T readScalarLine();

int main()
{

	std::cout << "Choose the scalar field of the components of the Matrix" << std::endl;
	std::cout << "(" << INT << " for ints, " << DOUBLE << " for double or " 
		  << COMPLEX << " for Complex):" << std::endl;
	getline(std::cin, g_line);
	int fieldChoice = atoi(g_line.c_str());
	assert (fieldChoice > 0 && fieldChoice <= NUM_KNOWN_FIELDS);

	std::cout << "You chose the field: " << FIELD_NAMES[fieldChoice-1] << std::endl;
	std::cout << "Choose operation:" << std::endl;

	int item;
	for (item = 1; item <= NUM_OPS; item++)
	{
		std::cout << item << ". " << OPERATIONS[item - 1] << std::endl;
	}


	getline(std::cin, g_line);
	int operatorChoice = atoi(g_line.c_str());
	assert (operatorChoice > 0 && operatorChoice <= NUM_OPS);

	if (operatorChoice > 2)
	{
		switch (fieldChoice)
		{
		 case INT:
			 unaryOperation<int>(operatorChoice);
			 break;
		 case DOUBLE:
			 unaryOperation<double>(operatorChoice);
			 break;
		 case COMPLEX:
			 unaryOperation<Complex>(operatorChoice);
			 break;
		}
	}
	else
	{
		switch (fieldChoice)
		{
		 case INT:
			 binaryOperation<int>(operatorChoice);
			 break;
		 case DOUBLE:
			 binaryOperation<double>(operatorChoice);
			 break;
		 case COMPLEX:
			 binaryOperation<Complex>(operatorChoice);
			 break;
		}
	}

	return 0;
}

template <typename T>
void checkIterators(Matrix<T>& mat){
std::cout<<"Vals:"<<std::endl;
	typename Matrix<T>::const_iterator it;
  for ( it = mat.begin() ; it != mat.end(); ++it){
    std::cout << ' ' << *it;
 }
  std::cout << '\n';

/* No need to imlement non-const  iterator 
  for (typename Matrix<T>::iterator it = mat.begin() ; it != mat.end(); ++it){
    (*it) = *it +1;;
 }
 std::cout <<std::endl << "Incremented Matrix:" <<std::endl <<mat <<std::endl;
 */
}

template <typename T>
void unaryOperation(const int operatorChoice)
{
	std::cout << "Operation " << OPERATIONS[operatorChoice-1] 
		  << " requires 1 operand Matrix." << std::endl;

	// Read the Matrix information:
	int rows, cols;
	std::vector<T> cells;
	readMatrixInfo(rows, cols, cells);
	Matrix<T> m(rows, cols, cells);
	std::cout << MAT_LINE << std::endl << "got Matrix:" << std::endl;
	std::cout << m;

	Matrix<T> resultMat;

	switch (operatorChoice)
	{
	 case TRANSPOSE:
		 resultMat = m.trans();
		 printResultMatrix(resultMat);
		 break;
/*	 case TRACE:
		 if (m.isSquareMatrix())
		 {
			 std::cout << "The Matrix is square and its trace is " << m.trace() << std::endl;
		 }
		 else
		 {
			 std::cout << "The Matrix isn't square" << std::endl;
		 }
		 break;
	 case SCALAR_TIMES_Matrix:
		 resultMat = scalar * m;
		 printResultMatrix(resultMat);
		 break;
	 case Matrix_TIMES_SCALAR:
		 resultMat = m * scalar;
		 printResultMatrix(resultMat);
		 break;
		 */
	}
	checkIterators(resultMat);
}


template <typename T>
void binaryOperation(const int operatorChoice)
{
	std::cout << "Operation " << OPERATIONS[operatorChoice-1] 
		  << " requires 2 operand matrices." << std::endl;

	// Read the Matrix information:
	int rows1, cols1, rows2, cols2;
	std::vector<T> cells1;
	std::vector<T> cells2;

	std::cout << "Insert first Matrix:" << std::endl;
	readMatrixInfo(rows1, cols1, cells1);
	std::cout << "Insert second Matrix:" << std::endl;
	readMatrixInfo(rows2, cols2, cells2);

	Matrix<T> m1(rows1, cols1, cells1);
	Matrix<T> m2(rows2, cols2, cells2);

	std::cout << MAT_LINE << std::endl << "Got first Matrix:" << std::endl;
	std::cout << m1;

	std::cout << MAT_LINE << std::endl << "Got second Matrix:" << std::endl;
	std::cout << m2;

	Matrix<T> resultMat;
assert(m1==m1);

	switch(operatorChoice)
	{
	 case ADD:
		 try
		 {
			 resultMat = m1 + m2;
			
		 }
		 catch (std::exception& exception)
		 {
			 std::cout << "Got Exception from Matrix with message: "
				   << std::endl << exception.what() << std::endl;
			 exit(1);
		 }
		 break;
	 case MUL:
		 try
		 {

			 assert(m1.cols() == m2.rows());
			 resultMat = m1 * m2;
		 }
		 catch (std::exception& exception)
		 {
			 std::cout << "Got Exception from Matrix with message: "
				   << std::endl << exception.what() << std::endl;
			 exit(1);
		 }
		 break;
	}
	printResultMatrix(resultMat);
	checkIterators(resultMat);
	}

template <typename T>
void readMatrixInfo(int& rows, int& cols, std::vector<T>& cells)
{
	std::cout << "number of myRows:";
	getline(std::cin, g_line);
	rows = atoi(g_line.c_str());

	std::cout << "number of columns:";
	getline(std::cin, g_line);
	cols = atoi(g_line.c_str());

	std::cout << "Now insert the values of the Matrix, row by row." << std::endl <<
		"After each cell add the char \'" << DELIM << "\'" << 
		" (including after the last cell of a row)." << std::endl << 
		"Each row should be in a separate line." << std::endl;

	int row, col;
	for (row = 0; row < rows; row++)
	{
		getline(std::cin, g_line);

		std::stringstream stream(g_line);
		std::string numStr;
		T val;
		for (col = 0; col < cols; col++)
		{
			getline(stream, numStr, DELIM);
			getNumFromString(numStr, &val);
			cells.push_back(val);
		}
	}

}

template <typename T>
void getNumFromString(const std::string &str, T *num)
{
	 //istringstream iss(s);
	T number(str);

	*num = number;
}




void getNumFromString(const std::string &str, Complex *num)
{
	 std::istringstream iss(str);
	 double r,i=0;
	 iss>>r;
//	 bool b=(iss>>i);
	 if (!iss ){
	 	i=0;
	 };
	Complex number(r,i);
	*num = number;
}

template <>
void getNumFromString(const std::string &str, int *num)
{
	int number = atoi(str.c_str());

	*num = number;
}
template <>
void getNumFromString(const std::string &str, double *num)
{
	double number = atof(str.c_str());

	*num = number;
}

template <typename T>
const T readScalarLine()
{
	getline(std::cin, g_line);
	T number;
	getNumFromString(g_line, &number);

	return number;
}


template <typename T>
void printResultMatrix(const Matrix<T>& mat)
{
	std::cout << LINE << std::endl;
	std::cout << "Resulted Matrix:" << std::endl;
	std::cout << mat;
}

