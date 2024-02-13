#include <vector>
#include <tuple>
#include <stdexcept>
#include <string>
#include <limits>
#include <type_traits>
#include <omp.h>

#pragma once

#define AXIS_X 0
#define AXIS_Y 1
#define ALL -1

#define HORIZONTAL 0
#define VERTICAL 1

// The axis about which certain operations will occur
typedef int axis_t;

// The memory orientation of the vector - used for SIMD optimisation.
typedef int m_orientation;

/**
 * This matrix class provides defines the mathematical behaviour of
 * a matrix.
 */
template<typename T>
class Matrix2D_S {

public:

	// This block prevents non-numerical matrix types.
	static_assert(std::is_arithmetic<T>::value, "Illegal type, template type T must be numerical.");

	/**
	 * Constructor - copies data from a vector
	 * 
	 * @param data: A 2D vector containg the input data
	 * 
	 * @param orientation: The memory orientation of the matrix being created.
	 */
	Matrix2D_S(std::vector<std::vector<T>> data, m_orientation orientation);

	/**
	 * Constructor - creates a matrix of a given size.
	 * 
	 * @param rows: The number of rows in the matrix being instantiated.
	 * 
	 * @param cols: The number of columsn int the matrix being instantiated.
	 * 
	 * @param orientation: The memory orientation of the matrix being created.
	 */
	Matrix2D_S(std::size_t rows, std::size_t cols, m_orientation orientation);

	/**
	 * Constructor - creates a matrix of a given size filled with a single value.
	 * 
	 * @param rows: The number of rows in the matrix being instantiated.
	 * 
	 * @param cols: The number of columsn int the matrix being instantiated.
	 * 
	 * @param orientation: The memory orientation of the matrix being created.
	 * 
	 * @param initializer: The value with which to initialize the matrix.
	 */
	Matrix2D_S(std::size_t rows, std::size_t cols, m_orientation orientation, T initializer);

	/**
	 * Constructor - creates a matrix of a given size filled using a generator function.
	 * 
	 * @param rows: The number of rows in the matrix being instantiated.
	 * 
	 * @param cols: The number of columsn int the matrix being instantiated.
	 * 
	 * @param orientation: The memory orientation of the matrix being created.
	 * 
	 * @param operation: The function with which to initialize the matrix.
	 */
	Matrix2D_S(std::size_t rows, std::size_t cols, m_orientation orientation, T(*operation)());

	/**
	 * Destructor
	 */
	~Matrix2D_S();

	/**
	 * This method finds the dot product of two matrices - output has horizontally
	 * aligned memory.
	 * 
	 * @param matrix: The matrix with which the dot product is being performed.
	 */
	Matrix2D_S<T> dot_H(Matrix2D_S<T>& matrix);

	/**
	* This method finds the dot product of two matrices - output has vertically
	* aligned memory.
	* 
	* @param matrix: The matrix with which the dot product is being performed.
	*/
	Matrix2D_S<T> dot_V(Matrix2D_S<T>& matrix);

	/**
	* This method finds the multiplication of a salar and a matrice	- output has
	* horizontally aligned memory.
	* 
	* @param scalar: The scalar with which to multiply this matrix.
	*/
	Matrix2D_S<T> multiply_H(T scalar);

	/**
	* This method finds the multiplication of a salar and a matrice	- output has
	* vertically aligned memory.
	* 
	* @param scalar: The scalar with which to multiply this matrix.
	*/
	Matrix2D_S<T> multiply_V(T scalar);

	/**
	* This method finds the multiplication of a salar and a matrice	- output has
	* horizontally aligned memory.
	* 
	* @param scalar: The scalar with which to divide this matrix.
	*/
	Matrix2D_S<T> divide_H(T scalar);

	/**
	* This method finds the multiplication of a salar and a matrice	- output has
	* vertically aligned memory.
	* 
	* @param scalar: The scalar with which to divide this matrix.
	*/
	Matrix2D_S<T> divide_V(T scalar);

	/**
	* This method finds the addition of a scalar to a matrice - output has
	* horizontally aligned memory.
	* 
	* @param scalar: The scalar with which to add this matrix.
	*/
	Matrix2D_S<T> add_H(T scalar);

	/**
	* This method finds the addition of a scalar to a matrice - output has
	* vertically aligned memory.
	* 
	* @param scalar: The scalar with which to add this matrix.
	*/
	Matrix2D_S<T> add_V(T scalar);

	/**
	* This method finds the addition of two matrices - output has
	* horizontally aligned memory.
	* 
	* @param matrix: The matrix with which to add this matrix.
	*/
	Matrix2D_S<T> add_H(Matrix2D_S<T>& matrix);

	/**
	* This method finds the addition of two matrices - output has
	* vertically aligned memory.
	* 
	* @param matrix: The matrix with which to add this matrix.
	*/
	Matrix2D_S<T> add_V(Matrix2D_S<T>& matrix);

	/**
	* This method finds the subtaction of a scalar to a matrice - output has
	* horizontally aligned memory.
	* 
	* @param scalar: The scalar with which to subract this matrix.
	*/
	Matrix2D_S<T> subtract_H(T scalar);

	/**
	* This method finds the subtaction of a scalar to a matrice - output has
	* vertically aligned memory.
	* 
	* @param scalar: The scalar with which to subract this matrix.
	*/
	Matrix2D_S<T> subtract_V(T scalar);

	/**
	* This method finds the subtraction of a matrice - output has horizontally
	* aligned memory.
	* 
	* @param matrix: The matrix with which to subtract this matrix.
	*/
	Matrix2D_S<T> subtract_H(Matrix2D_S<T>& matrix);

	/**
	* This method finds the subtraction of a matrice - output has vertically
	* aligned memory.
	* 
	* @param matrix: The matrix with which to subtract this matrix.
	*/
	Matrix2D_S<T> subtract_V(Matrix2D_S<T>& matrix);

	/**
	* This method finds the elementwise max of a scalar - output has horizontally
	* aligned memory.
	*
	* @param scalar: The scalar against which the maximum is found.
	*/
	Matrix2D_S<T> max_H(T scalar);

	/**
	* This method finds the elementwise max of a scalar - output has vertically
	* aligned memory.
	* 
	* @param scalar: The scalar against which the maximum is found.
	*/
	Matrix2D_S<T> max_V(T scalar);

	/**
	* This method finds the elementwise min of a scalar - output has horizontally
	* aligned memory.
	* 
	* @param scalar: The scalar against which the minimum is found.
	*/
	Matrix2D_S<T> min_H(T scalar);

	/**
	* This method finds the elementwise min of a scalar - output has vertically
	* aligned memory.
	* 
	* @param scalar: The scalar against which the minimum is found.
	*/
	Matrix2D_S<T> min_V(T scalar);

	/**
	* This method applies a single operation to the entire matrix elementwise.
	* The output is stored in horizontally aligned memory.
	* 
	* @param operation: The operation which is performed on the entire array.
	*/
	Matrix2D_S<T> elementwiseApply_H(T(*operation)(T));

	/**
	* This method applies a single operation to the entire matrix elementwise.
	* The output is stored in vertically aligned memory.
	*
	* @param operation: The operation which is performed on the entire array.
	*/
	Matrix2D_S<T> elementwiseApply_V(T(*operation)(T));

	/**
	* This method returns the Matrix2D_S transform of the given matrice.
	*/
	void transform();

	/**
	* This method re-allocates the memory so that it is aligned perpendicularly
	* to how it was previously.
	*/
	void realign();

	/**
	 * This getter method returns a tuple of the array's dimensions.
	 */
	std::tuple<std::size_t, std::size_t> dimensions();

	/**
	 * This getter method returns the memory orientation of the array.
	 */
	m_orientation getOrientation();

	/**
	 * This getter method returns a reference to the matrix data (vector).
	 */
	T** getMatrix();

	/**
	 * This method converts the Matrix to a string.
	 */
	std::string toString();

private:

	/**
	 * The number of rows.
	 */
	std::size_t rows;

	/**
	 * The number of columns.
	 */
	std::size_t cols;

	/**
	 * The total number of elements in the matrix.
	 */
	std::size_t size;

	/**
	 * The orientation of memory as assigned for the Matrix.
	 */
	m_orientation orientation;

	/**
	 * This nested array holds the matrix data.
	 */
	T** matrix;

	/**
	 * This static value holds the maximum value possible for type T
	 */
	static T maximum;

	/**
	 * This static value holds the minimum value possible for type T
	 */
	static T minimum;

	/**
	 * This method performs creates a new empty matrix - output has horizontally
	 * aligned memory.
	 * 
	 * @param rows: The number of rows in the matrix being instantiated.
	 * 
	 * @param cols: The number of columsn int the matrix being instantiated.
	 * 
	 * @param size: The total number of elements in the array.
	 */
	T** create_H(std::size_t rows, std::size_t cols, std::size_t size);

	/**
	* This method performs creates a new empty matrix - output has vertically
	* aligned memory.
	* 
	* @param rows: The number of rows in the matrix being instantiated.
	* 
	* @param cols: The number of columsn int the matrix being instantiated.
	* 
	* @param size: The total number of elements in the array.
	*/
	T** create_V(std::size_t rows, std::size_t cols, std::size_t size);

	/**
	* Constructor - creates a matrix from a pointer - dangerous private method.
	* 
	* @param data: The pointer to the array of data.
	* 
	* @param rows: The number of rows in the matrix being instantiated.
	* 
	* @param cols: The number of columsn int the matrix being instantiated.
	* 
	* @param size: The total number of elements in the array.
	* 
	* @param orientation: The memory orientation of the matrix being created.
	*/
	Matrix2D_S(T** data, std::size_t rows, std::size_t cols, std::size_t size, m_orientation orientation);
};




template<typename T>
T Matrix2D_S<T>::minimum = std::numeric_limits<T>::min();




template<typename T>
T Matrix2D_S<T>::maximum = std::numeric_limits<T>::max();




template<typename T>
T** Matrix2D_S<T>::create_H(std::size_t rows, std::size_t cols, std::size_t size) {

	// Throw an exception for a matrix size of 0.
	if (size == 0) {
		throw std::invalid_argument("Invalid matrix size, the total size (rows x columns) must be greater than 0.");
	}

	// The row pointers are created
	T** matrix = new T * [rows];

	// All the data is allocated contiguously.
	matrix[0] = new T[size];

	// Pointers are set to the contiguous block.
	for (std::size_t i = 1; i < rows; i++) {
		matrix[i] = &matrix[0][cols * i];
	}

	return matrix;
}




template<typename T>
T** Matrix2D_S<T>::create_V(std::size_t rows, std::size_t cols, std::size_t size) {

	// Throw an exception for a matrix size of 0.
	if (size == 0) {
		throw std::invalid_argument("Invalid matrix size, the total size (rows x columns) must be greater than 0.");
	}

	// The column pointers are created
	 T** matrix = new T * [cols];

	// All the data is allocated contiguously.
	matrix[0] = new T[size];

	// Pointers are set to the contiguous block.
	for (std::size_t i = 1; i < cols; i++) {
		matrix[i] = &matrix[0][rows * i];
	}

	return matrix;
}




template<typename T>
Matrix2D_S<T>::Matrix2D_S(std::vector<std::vector<T>> data, m_orientation orientation) {

	this->rows = data.size();
	this->cols = data.at(0).size();
	this->orientation = orientation;
	this->size = this->rows * this->cols;

	// The memory is allocated to the heap.
	if (orientation == HORIZONTAL) {

		this->matrix = this->create_H(this->rows, this->cols, this->size);

		// Data is copied from a vector to an array.
		for (std::size_t i = 0; i < this->rows; i++) {

			if (this->cols != data.at(i).size()) {
				throw std::invalid_argument("Dimension inconsistencies were found in the input data.");
			}

			for (std::size_t j = 0; j < this->cols; j++) {

				this->matrix[i][j] = data[i][j];

			}
		}

	} else if(orientation == VERTICAL) {

		this->matrix = this->create_V(this->rows, this->cols, this->size);

		// Data is copied from a vector to an array.
		for (std::size_t i = 0; i < this->rows; i++) {

			if (this->cols != data.at(i).size()) {
				throw std::invalid_argument("Dimension inconsistencies were found in the input data.");
			}

			for (std::size_t j = 0; j < this->cols; j++) {

				this->matrix[j][i] = data[i][j];

			}
		}

	} else {
		throw std::invalid_argument("No such orientation exists.");
	}
}




template<typename T>
Matrix2D_S<T>::Matrix2D_S(std::size_t rows, std::size_t cols, m_orientation orientation) {

	this->rows = rows;
	this->cols = cols;
	this->orientation = orientation;
	this->size = this->rows * this->cols;

	// The memory is allocated to the heap.
	if (orientation == HORIZONTAL) {
		this->matrix = this->create_H(this->rows, this->cols, this->size);
	} else if (orientation == VERTICAL) {
		this->matrix = this->create_V(this->rows, this->cols, this->size);
	} else {
		throw std::invalid_argument("No such orientation exists.");
	}
}




template<typename T>
Matrix2D_S<T>::Matrix2D_S(std::size_t rows, std::size_t cols, m_orientation orientation, T initializer) {

	this->rows = rows;
	this->cols = cols;
	this->orientation = orientation;
	this->size = this->rows * this->cols;

	// The memory is allocated to the heap.
	if (orientation == HORIZONTAL) {
		this->matrix = this->create_H(this->rows, this->cols, this->size);
	} else if (orientation == VERTICAL) {
		this->matrix = this->create_V(this->rows, this->cols, this->size);
	} else {
		throw std::invalid_argument("No such orientation exists.");
	}

	// Matrix is populated
	for (std::size_t i = 0; i < this->size; i++) {
		this->matrix[0][i] = initializer;
	}
}




template<typename T>
Matrix2D_S<T>::Matrix2D_S(std::size_t rows, std::size_t cols, m_orientation orientation, T(*operation)()) {

	this->rows = rows;
	this->cols = cols;
	this->orientation = orientation;
	this->size = this->rows * this->cols;

	// The memory is allocated to the heap.
	if (orientation == HORIZONTAL) {
		this->matrix = this->create_H(this->rows, this->cols, this->size);
	} else if (orientation == VERTICAL) {
		this->matrix = this->create_V(this->rows, this->cols, this->size);
	} else {
		throw std::invalid_argument("No such orientation exists.");
	}

	// Matrix is populated
	for (std::size_t i = 0; i < this->size; i++) {
		this->matrix[0][i] = operation();
	}
}



template<typename T>
Matrix2D_S<T>::Matrix2D_S(T** data, std::size_t rows, std::size_t cols, std::size_t size, m_orientation orientation) {

	this->rows = rows;
	this->cols = cols;
	this->orientation = orientation;
	this->size = this->rows * this->cols;
	this->matrix = data;
}




template<typename T>
Matrix2D_S<T>::~Matrix2D_S() {
	// Delete old array.
	delete[] this->matrix[0];
	delete[] this->matrix;
}




template<typename T>
Matrix2D_S<T> Matrix2D_S<T>::dot_H(Matrix2D_S<T>& matrix) {

	std::size_t numRows = std::get<0>(matrix.dimensions());
	std::size_t numCols = std::get<1>(matrix.dimensions());

	// Conditional throws exception for incompatible dimensions or orientations
	if (this->cols != numRows) {
		throw std::invalid_argument("The number of columns in matrix 1 does not equal the number of rows in matrix 2.");
	} else if (this->getOrientation() != HORIZONTAL || matrix.getOrientation() != VERTICAL) {
		throw std::invalid_argument("The orientation of matrix 1 must be horizontal, and the orientation of matrix 2 "
									"must be vertical - for simd optimization.");
	}

	T** externalMatrix = matrix.getMatrix();

	std::size_t size = this->rows * numCols;

	// The output matrix is created
	T** outputMatrix = this->create_H(this->rows, numCols, size);
	Matrix2D_S<T> output(outputMatrix, this->rows, numCols, size, HORIZONTAL);

	T sum;

	// For each row (internal)
	for (std::size_t i = 0; i < this->rows; i++) {
		// Multiply and sum by all columns (external)
		for (std::size_t j = 0; j < numCols; j++) {

			sum = 0;			

			// For each row (external)
			#pragma omp simd reduction(+:sum)
			for (std::size_t k = 0; k < numRows; k++) {

				sum += this->matrix[i][k] * externalMatrix[j][k];

			}

			outputMatrix[i][j] = sum;
		}
	}

	return output;
}




template<typename T>
Matrix2D_S<T> Matrix2D_S<T>::dot_V(Matrix2D_S<T>& matrix) {

	std::size_t numRows = std::get<0>(matrix.dimensions());
	std::size_t numCols = std::get<1>(matrix.dimensions());

	// Conditional throws exception for incompatible dimensions or orientations
	if (this->cols != numRows) {
		throw std::invalid_argument("The number of columns in matrix 1 does not equal the number of rows in matrix 2.");
	} else if (this->getOrientation() != HORIZONTAL || matrix.getOrientation() != VERTICAL) {
		throw std::invalid_argument("The orientation of matrix 1 must be horizontal, and the orientation of matrix 2 "
									"must be vertical - for simd optimization.");
	}

	T** externalMatrix = matrix.getMatrix();

	std::size_t size = this->rows * numCols;

	// The output matrix is created
	T** outputMatrix = this->create_V(this->rows, numCols, size);
	Matrix2D_S<T> output(outputMatrix, this->rows, numCols, size, VERTICAL);

	T sum = 0;

	// For each row (internal)
	for (std::size_t i = 0; i < this->rows; i++) {
		// Multiply and sum by all columns (external)
		for (std::size_t j = 0; j < numCols; j++) {

			sum = 0;

			// For each row (external)
			#pragma omp simd reduction(+:sum)
			for (std::size_t k = 0; k < numRows; k++) {

				sum += this->matrix[i][k] * externalMatrix[j][k];

			}

			outputMatrix[j][i] = sum;
		}
	}

	return output;
}





template<typename T>
Matrix2D_S<T> Matrix2D_S<T>::multiply_H(T scalar) {

	// The output matrix is created
	T** outputMatrix = this->create_H(this->rows, this->cols, this->size);
	Matrix2D_S<T> output(outputMatrix, this->rows, this->cols, this->size, HORIZONTAL);

	// Multiply each element by the scalar
	for (std::size_t i = 0; i < size; ++i) {
		outputMatrix[0][i] = this->matrix[0][i] * scalar;
	}

	return output;
}




template<typename T>
Matrix2D_S<T> Matrix2D_S<T>::multiply_V(T scalar) {

	// The output matrix is created
	T** outputMatrix = this->create_V(this->rows, this->cols, this->size);
	Matrix2D_S<T> output(outputMatrix, this->rows, this->cols, this->size, VERTICAL);

	// Multiply each element by the scalar
	for (std::size_t i = 0; i < size; ++i) {
		outputMatrix[0][i] = this->matrix[0][i] * scalar;
	}

	return output;
}




template<typename T>
Matrix2D_S<T> Matrix2D_S<T>::divide_H(T scalar) {

	// The output matrix is created
	T** outputMatrix = this->create_H(this->rows, this->cols, this->size);
	Matrix2D_S<T> output(outputMatrix, this->rows, this->cols, this->size, HORIZONTAL);

	// Data is divided
	for (std::size_t i = 0; i < size; i++) {
		outputMatrix[0][i] = this->matrix[0][i] / scalar;
	}

	return output;
}




template<typename T>
Matrix2D_S<T> Matrix2D_S<T>::divide_V(T scalar) {

	// The output matrix is created
	T** outputMatrix = this->create_V(this->rows, this->cols, this->size);
	Matrix2D_S<T> output(outputMatrix, this->rows, this->cols, this->size, VERTICAL);

	// Data is divided
	for (std::size_t i = 0; i < size; i++) {
		outputMatrix[0][i] = this->matrix[0][i] / scalar;
	}

	return output;
}




template<typename T>
Matrix2D_S<T> Matrix2D_S<T>::add_H(T scalar) {

	// The output matrix is created
	T** outputMatrix = this->create_H(this->rows, this->cols, this->size);
	Matrix2D_S<T> output(outputMatrix, this->rows, this->cols, this->size, HORIZONTAL);

	// Data is added
	for (std::size_t i = 0; i < size; i++) {
		outputMatrix[0][i] = this->matrix[0][i] + scalar;
	}

	return output;
}




template<typename T>
Matrix2D_S<T> Matrix2D_S<T>::add_V(T scalar) {

	// The output matrix is created
	T** outputMatrix = this->create_V(this->rows, this->cols, this->size);
	Matrix2D_S<T> output(outputMatrix, this->rows, this->cols, this->size, VERTICAL);

	// Data is added
	for (std::size_t i = 0; i < size; i++) {
		outputMatrix[0][i] = this->matrix[0][i] + scalar;
	}

	return output;
}




template<typename T>
Matrix2D_S<T> Matrix2D_S<T>::add_H(Matrix2D_S<T>& matrix) {

	std::size_t numRows = std::get<0>(matrix.dimensions());
	std::size_t numCols = std::get<1>(matrix.dimensions());

	if (this->rows != numRows || this->cols != numCols) {
		throw std::invalid_argument("These matrices cannot be added, the dimensions are unequal.");
	} else if (this->getOrientation() != matrix.getOrientation()) {
		throw std::invalid_argument("The orientation of the two matrices must be the same.");
	}

	T** inputMatrix = matrix.getMatrix();

	// The output matrix is created
	T** outputMatrix = this->create_H(this->rows, this->cols, this->size);
	Matrix2D_S<T> output(outputMatrix, this->rows, this->cols, this->size, HORIZONTAL);

	// Data is added
	for (std::size_t i = 0; i < this->size; i++) {
		outputMatrix[0][i] = this->matrix[0][i] + inputMatrix[0][i];
	}

	return output;
}




template<typename T>
Matrix2D_S<T> Matrix2D_S<T>::add_V(Matrix2D_S<T>& matrix) {

	std::size_t numRows = std::get<0>(matrix.dimensions());
	std::size_t numCols = std::get<1>(matrix.dimensions());

	if (this->rows != numRows || this->cols != numCols) {
		throw std::invalid_argument("These matrices cannot be added, the dimensions are unequal.");
	} else if (this->getOrientation() != matrix.getOrientation()) {
		throw std::invalid_argument("The orientation of the two matrices must be the same.");
	}

	T** inputMatrix = matrix.getMatrix();

	// The output matrix is created
	T** outputMatrix = this->create_V(this->rows, this->cols, this->size);
	Matrix2D_S<T> output(outputMatrix, this->rows, this->cols, this->size, VERTICAL);

	// Data is added
	for (std::size_t i = 0; i < size; i++) {
		outputMatrix[0][i] = this->matrix[0][i] + inputMatrix[0][i];
	}

	return output;
}




template<typename T>
Matrix2D_S<T> Matrix2D_S<T>::subtract_H(T scalar) {

	// The output matrix is created
	T** outputMatrix = this->create_H(this->rows, this->cols, this->size);
	Matrix2D_S<T> output(outputMatrix, this->rows, this->cols, this->size, HORIZONTAL);

	// Data is subtracted
	for (std::size_t i = 0; i < size; i++) {
		outputMatrix[0][i] = this->matrix[0][i] - scalar;
	}

	return output;
}





template<typename T>
Matrix2D_S<T> Matrix2D_S<T>::subtract_V(T scalar) {

	// The output matrix is created
	T** outputMatrix = this->create_V(this->rows, this->cols, this->size);
	Matrix2D_S<T> output(outputMatrix, this->rows, this->cols, this->size, VERTICAL);

	// Data is subtracted
	for (std::size_t i = 0; i < this->size; i++) {
		outputMatrix[0][i] = this->matrix[0][i] - scalar;
	}

	return output;
}




template<typename T>
Matrix2D_S<T> Matrix2D_S<T>::subtract_H(Matrix2D_S<T>& matrix) {

	std::size_t numRows = std::get<0>(matrix.dimensions());
	std::size_t numCols = std::get<1>(matrix.dimensions());

	if (this->rows != numRows || this->cols != numCols) {
		throw std::invalid_argument("These matrices cannot be subtracted, the dimensions are unequal.");
	} else if (this->getOrientation() != matrix.getOrientation()) {
		throw std::invalid_argument("The orientation of the two matrices must be the same.");
	}

	T** inputMatrix = matrix.getMatrix();

	// Data is subtracted
	T** outputMatrix = this->create_H(this->rows, this->cols, this->size);
	Matrix2D_S<T> output(outputMatrix, this->rows, this->cols, this->size, HORIZONTAL);

	// Data is added
	for (std::size_t i = 0; i < this->size; i++) {
		outputMatrix[0][i] = this->matrix[0][i] - inputMatrix[0][i];
	}

	return output;
}





template<typename T>
Matrix2D_S<T> Matrix2D_S<T>::subtract_V(Matrix2D_S<T>& matrix) {

	std::size_t numRows = std::get<0>(matrix.dimensions());
	std::size_t numCols = std::get<1>(matrix.dimensions());

	if (this->rows != numRows || this->cols != numCols) {
		throw std::invalid_argument("These matrices cannot be subtracted, the dimensions are unequal.");
	} else if (this->getOrientation() != matrix.getOrientation()) {
		throw std::invalid_argument("The orientation of the two matrices must be the same.");
	}

	// The output matrix is created
	T** outputMatrix = this->create_V(this->rows, this->cols, this->size);
	Matrix2D_S<T> output(outputMatrix, this->rows, this->cols, this->size, VERTICAL);

	// Data is subtracted
	for (std::size_t i = 0; i < this->size; i++) {
		outputMatrix[0][i] = this->matrix[0][i] - inputMatrix[0][i];
	}

	return output;
}





template<typename T>
Matrix2D_S<T> Matrix2D_S<T>::min_H(T scalar) {

	// The output matrix is created
	T** outputMatrix = this->create_H(this->rows, this->cols, this->size);
	Matrix2D_S<T> output(outputMatrix, this->rows, this->cols, this->size, HORIZONTAL);

	// Elementwise maximum is found
	for (std::size_t i = 0; i < size; i++) {
		outputMatrix[0][i] = (this->matrix[0][i] < scalar ? this->matrix[0][i] : scalar);
	}

	return output;
}






template<typename T>
Matrix2D_S<T> Matrix2D_S<T>::min_V(T scalar) {

	// The output matrix is created
	T** outputMatrix = this->create_V(this->rows, this->cols, this->size);
	Matrix2D_S<T> output(outputMatrix, this->rows, this->cols, this->size, VERTICAL);

	// Elementwise maximum is found
	for (std::size_t i = 0; i < size; i++) {
		outputMatrix[0][i] = (this->matrix[0][i] < scalar ? this->matrix[0][i] : scalar);
	}

	return output;
}






template<typename T>
Matrix2D_S<T> Matrix2D_S<T>::max_H(T scalar) {

	// The output matrix is created
	T** outputMatrix = this->create_H(this->rows, this->cols, this->size);
	Matrix2D_S<T> output(outputMatrix, this->rows, this->cols, this->size, HORIZONTAL);

	// Elementwise maximum is found
	for (std::size_t i = 0; i < size; i++) {
		outputMatrix[0][i] = (this->matrix[0][i] > scalar ? this->matrix[0][i] : scalar);
	}

	return output;
}






template<typename T>
Matrix2D_S<T> Matrix2D_S<T>::max_V(T scalar) {

	// The output matrix is created
	T** outputMatrix = this->create_V(this->rows, this->cols, this->size);
	Matrix2D_S<T> output(outputMatrix, this->rows, this->cols, this->size, VERTICAL);

	// Elementwise maximum is found
	for (std::size_t i = 0; i < size; i++) {
		outputMatrix[0][i] = (this->matrix[0][i] > scalar ? this->matrix[0][i] : scalar);
	}

	return output;
}






template<typename T>
Matrix2D_S<T> Matrix2D_S<T>::elementwiseApply_H(T(*operation)(T)) {

	// The output matrix is created
	T** outputMatrix = this->create_H(this->rows, this->cols, this->size);
	Matrix2D_S<T> output(outputMatrix, this->rows, this->cols, this->size, HORIZONTAL);

	// Elementwise maximum is found
	for (std::size_t i = 0; i < size; i++) {
		outputMatrix[0][i] = operation(this->matrix[0][i]);
	}

	return output;
}





template<typename T>
Matrix2D_S<T> Matrix2D_S<T>::elementwiseApply_V(T(*operation)(T)) {

	// The output matrix is created
	T** outputMatrix = this->create_V(this->rows, this->cols, this->size);
	Matrix2D_S<T> output(outputMatrix, this->rows, this->cols, this->size, VERTICAL);

	// Elementwise maximum is found
	for (std::size_t i = 0; i < size; i++) {
		outputMatrix[0][i] = operation(this->matrix[0][i]);
	}

	return output;
}





template<typename T>
void Matrix2D_S<T>::transform() {

	T** replaceMatrix = nullptr;

	if (this->orientation == HORIZONTAL) {
		replaceMatrix = this->create_H(this->cols, this->rows, this->size);
	} else {
		replaceMatrix = this->create_V(this->cols, this->rows, this->size);
	}

	// Data is transformed (rows,cols) => (cols, rows)
	for (std::size_t i = 0; i < this->cols; i++) {
		#pragma omp simd
		for (std::size_t j = 0; j < this->rows; j++) {

			replaceMatrix[i][j] = this->matrix[j][i];

		}
	}

	// Delete old array.
	delete[] this->matrix[0];
	delete[] this->matrix;

	// Assign transformed array.
	this->matrix = replaceMatrix;

	// Set the row and column values
	std::size_t temp = this->rows;
	this->rows = this->cols;
	this->cols = temp;
}




template<typename T>
void Matrix2D_S<T>::realign() {

	T** replaceMatrix = nullptr;

	if (this->orientation == HORIZONTAL) {
		replaceMatrix = this->create_V(this->rows, this->cols, this->size);
		this->orientation = VERTICAL;
	} else {
		replaceMatrix = this->create_H(this->rows, this->cols, this->size);
		this->orientation = HORIZONTAL;
	}

	// Data is transferred
	for (std::size_t i = 0; i < this->size; i++) {
		replaceMatrix[0][i] = this->matrix[0][i];
	}

	// Delete old array.
	delete[] this->matrix[0];
	delete[] this->matrix;

	// Assign transformed array.
	this->matrix = replaceMatrix;	
}




template<typename T>
std::tuple<std::size_t, std::size_t> Matrix2D_S<T>::dimensions() {

	std::tuple<std::size_t, std::size_t> dim{ this->rows, this->cols };

	return dim;
}





template<typename T>
m_orientation Matrix2D_S<T>::getOrientation() {
	return this->orientation;
}





template<typename T>
T** Matrix2D_S<T>::getMatrix() {
	return this->matrix;
}





template<typename T>
std::string Matrix2D_S<T>::toString() {

	std::string matrix;

	if(orientation == HORIZONTAL){

		// Iterate over matrix and convert to string data.
		for (std::size_t i = 0; i < this->rows; i++) {
			for (std::size_t j = 0; j < this->cols; j++) {
				matrix += std::to_string(this->matrix[i][j]) + ", ";
			}
			matrix += "\n";
		}

	}else{
	
		// Iterate over matrix and convert to string data.
		for (std::size_t i = 0; i < this->rows; i++) {
			for (std::size_t j = 0; j < this->cols; j++) {
				matrix += std::to_string(this->matrix[j][i]) + ", ";
			}
			matrix += "\n";
		}

	}

	return matrix;
}