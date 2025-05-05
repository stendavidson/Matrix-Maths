#pragma once

#include <vector>
#include <stdexcept>
#include <string>
#include <limits>
#include <type_traits>
#include <immintrin.h>



// Check if AVX2 is available
#if (defined(__AVX__) && defined(__AVX2__)) || defined(_M_AVX2)
#define SIMD_SUPPORTED true
#else
#define SIMD_SUPPORTED false
#endif



// Static assert to ensure the required SIMD is available at compile time
static_assert(SIMD_SUPPORTED, "Required SIMD support (AVX2) is not available on this platform!");



// The memory allocation method is determined by the compiler
#if defined(_MSC_VER)
#include <malloc.h>
#include <intrin.h>

/**
 * This function will allocate contiguous memory with the sepcified alignment.
 *
 * @param alignment The alignment of the allocation.
 *
 * @param size The size of the allocation.
 *
 * @return A pointer to the newly allocated memory.
 */
inline void* aligned_alloc(size_t alignment, size_t size) {
	return _aligned_malloc(size, alignment);
}

/**
 * Free the memory at given pointer
 *
 * @param ptr The pointer to the memory to be de-allocated.
 */
inline void aligned_free(void* ptr) {
	_aligned_free(ptr);
}

#else
#include <cstdlib>
#include <cpuid.h>

/**
 * This function will allocate contiguous memory with the sepcified alignment.
 *
 * @param alignment The alignment of the allocation.
 *
 * @param size The size of the allocation.
 *
 * @return A pointer to the newly allocated memory.
 */
inline void* aligned_alloc(size_t alignment, size_t size) {
	return std::aligned_alloc(alignment, size);
}

/**
 * Free the memory at given pointer
 *
 * @param ptr The pointer to the memory to be de-allocated.
 */
inline void aligned_free(void* ptr) {
	std::free(ptr);
}

#endif

// Define axis of operations
#define AXIS_X 0
#define AXIS_Y 1
#define ALL -1

// Define alignments
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

	// This block prevents non-numerical and un-optimizable template data types
	static_assert(
		std::is_same<T, float>::value ||
		std::is_same<T, double>::value,
		"Template type must be one of: float and double."
		);

	// This block sets a secondary type for the intrinsics
	using intrinsic_type = typename std::conditional<std::is_same<T, float>::value, __m256, __m256d>::type;

	// This block sets a third type for masks
	using mask_type = typename std::conditional<std::is_same<T, float>::value, int, std::int64_t>::type;

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
	Matrix2D_S(long rows, long cols, m_orientation orientation);

	/**
	 * Copy Constructor
	 *
	 * @param matrix: The matrix being copied
	 */
	Matrix2D_S(const Matrix2D_S& matrix);

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
	Matrix2D_S(long rows, long cols, m_orientation orientation, T initializer);

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
	Matrix2D_S(long rows, long cols, m_orientation orientation, T(*operation)());

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
	* This method finds the multiplication of a salar and a matrice
	*
	* @param scalar: The scalar with which to multiply this matrix.
	*/
	Matrix2D_S<T> multiply(T scalar);

	/**
	* This method finds the elementwise multiplication of two matrices.
	*
	* @param scalar: The matrix with which to multiply this matrix.
	*/
	Matrix2D_S<T> multiply(Matrix2D_S<T>& matrix);

	/**
	* This method finds the multiplication of a salar and a matrice
	*
	* @param scalar: The scalar with which to divide this matrix.
	*/
	Matrix2D_S<T> divide(T scalar);

	/**
	* This method finds the elementwise division of two matrices.
	*
	* @param scalar: The matrix with which to divide this matrix.
	*/
	Matrix2D_S<T> divide(Matrix2D_S<T>& matrix);

	/**
	* This method finds the addition of a scalar to a matrice
	*
	* @param scalar: The scalar with which to add this matrix.
	*/
	Matrix2D_S<T> add(T scalar);

	/**
	* This method finds the addition of two matrices
	*
	* @param matrix: The matrix with which to add this matrix.
	*/
	Matrix2D_S<T> add(Matrix2D_S<T>& matrix);

	/**
	* This method finds the subtaction of a scalar to a matrice
	*
	* @param scalar: The scalar with which to subract this matrix.
	*/
	Matrix2D_S<T> subtract(T scalar);

	/**
	* This method finds the subtraction of a matrice.
	*
	* @param matrix: The matrix with which to subtract this matrix.
	*/
	Matrix2D_S<T> subtract(Matrix2D_S<T>& matrix);

	/**
	* This method finds the elementwise max of a scalar.
	*
	* @param scalar: The scalar against which the maximum is found.
	*/
	Matrix2D_S<T> max(T scalar);

	/**
	* This method finds the maximum values of two matrices.
	*
	* @param matrix: The matrix against which to find the maximum
	*/
	Matrix2D_S<T> max(Matrix2D_S<T>& matrix);

	/**
	* This method finds the elementwise min of a scalar
	*
	* @param scalar: The scalar against which the minimum is found.
	*/
	Matrix2D_S<T> min(T scalar);

	/**
	* This method finds the minimum values of two matrices.
	*
	* @param matrix: The matrix against which to find the minimum
	*/
	Matrix2D_S<T> min(Matrix2D_S<T>& matrix);

	/**
	* This method applies a single operation to the entire matrix elementwise.
	*
	* @param operation: The operation which is performed on the entire array.
	*/
	Matrix2D_S<T> apply(T(*operation)(T));

	/**
	* This method returns the Matrix2D_S transpose of the given matrice.
	*/
	Matrix2D_S<T> transpose();

	/**
	* This method re-allocates the memory so that it is aligned perpendicularly
	* to how it was previously.
	*/
	Matrix2D_S<T> realign();

	/**
	 * This getter method returns the the number of rows in the matrix
	 */
	long getRows();

	/**
	 * This getter method returns the the number of rows in the matrix
	 */
	long getCols();

	/**
	 * This getter method returns the amount of padding per row or column (depending on memory orientation)
	 */
	long getPadding();

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
	long rows;

	/**
	 * The number of columns.
	 */
	long cols;

	/**
	 * The total size of the array including padding
	 */
	long size;

	/**
	 * The column-wise or row-wise padding
	 */
	long padding;

	/**
	 * This static value holds the SIMD step size.
	 */
	static long simdStep;

	/**
	 * This static function holds the SIMD addition function
	 */
	static inline intrinsic_type intrinsic_add(intrinsic_type a, intrinsic_type b);

	/**
	 * This static function holds the SIMD aggregation funnction
	 */
	static inline T intrinsic_aggregate_add(intrinsic_type a);

	/**
	 * This static function holds the SIMD subtract function
	 */
	static inline intrinsic_type intrinsic_subtract(intrinsic_type a, intrinsic_type b);

	/**
	 * This static function holds the SIMD multiply function
	 */
	static inline intrinsic_type intrinsic_multiply(intrinsic_type a, intrinsic_type b);

	/**
	 * This static function holds the SIMD multiply AND add function
	 */
	static inline intrinsic_type intrinsic_multiply_add(intrinsic_type a, intrinsic_type b, intrinsic_type c);

	/**
	 * This static function holds the SIMD multiply function
	 */
	static inline intrinsic_type intrinsic_divide(intrinsic_type a, intrinsic_type b);

	/**
	 * This static function holds the SIMD elementwise minimum function
	 */
	static inline intrinsic_type intrinsic_minimum(intrinsic_type a, intrinsic_type b);

	/**
	 * This static function holds the SIMD elementwise maximum function
	 */
	static inline intrinsic_type intrinsic_maximum(intrinsic_type a, intrinsic_type b);

	/**
	 * This static function holds the SIMD elementwise bitwise AND function
	 */
	static inline intrinsic_type intrinsic_and(intrinsic_type a, intrinsic_type b);

	/**
	 * This static function holds the SIMD load function
	 */
	static inline intrinsic_type intrinsic_load(T* ptr);

	/**
	 * This static function holds a SIMD mask load function
	 */
	static inline intrinsic_type intrinsic_mask_load(mask_type* ptr);

	/**
	 * This static function holds the SIMD store function
	 */
	static inline void intrinsic_store(T* ptr, intrinsic_type a);

	/**
	 * This static function wraps the SIMD function that creats a zero initialized register
	 */
	static inline intrinsic_type intrinsic_zeroes_register();

	/**
	 * This static function wraps the SIMD function that creats register with all values set
	 * to a specific value.
	 */
	static inline intrinsic_type intrinsic_init_register(T a);

	/**
	 * The orientation of memory as assigned for the Matrix.
	 */
	m_orientation orientation;

	/**
	 * This nested array holds the matrix data.
	 */
	T** matrix;

	/**
	 * This array holds the row or column mask
	 */
	mask_type* mask;

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
	 *
	 * @param padding: The additional padding elements in the array.
	 */
	T** create_H(long rows, long cols, long size, long padding);

	/**
	* This method performs creates a new empty matrix - output has vertically
	* aligned memory.
	*
	* @param rows: The number of rows in the matrix being instantiated.
	*
	* @param cols: The number of columsn int the matrix being instantiated.
	*
	* @param size: The total number of elements in the array.
	*
	* @param padding: The additional padding elements in the array.
	*/
	T** create_V(long rows, long cols, long size, long padding);
};



template<typename T>
T Matrix2D_S<T>::minimum = std::numeric_limits<T>::min();



template<typename T>
T Matrix2D_S<T>::maximum = std::numeric_limits<T>::max();



template<>
long Matrix2D_S<float>::simdStep = 8;
template<>
long Matrix2D_S<double>::simdStep = 4;



template<>
inline typename Matrix2D_S<float>::intrinsic_type Matrix2D_S<float>::intrinsic_add(intrinsic_type a, intrinsic_type b) {
	return _mm256_add_ps(a, b);
}
template<>
inline typename Matrix2D_S<double>::intrinsic_type Matrix2D_S<double>::intrinsic_add(intrinsic_type a, intrinsic_type b) {
	return _mm256_add_pd(a, b);
}



template<>
inline float Matrix2D_S<float>::intrinsic_aggregate_add(intrinsic_type a) {

	// Array aggregated horizontally - for 8 values this is sees a performance improvement
	__m128 half1 = _mm256_castps256_ps128(a);
	__m128 half2 = _mm256_extractf128_ps(a, 1);
	__m128 sum128 = _mm_add_ps(half1, half2);

	// Now horizontally add the 4 floats in sum128.
	sum128 = _mm_hadd_ps(sum128, sum128);
	sum128 = _mm_hadd_ps(sum128, sum128);

	// Return the sum as a scalar
	return _mm_cvtss_f32(sum128);
}
template<>
inline double Matrix2D_S<double>::intrinsic_aggregate_add(intrinsic_type a) {

	#if defined(_MSC_VER)
	__declspec(align(32)) double data[4];
	#else
	double data[4] __attribute__((aligned(32)));
	#endif

	_mm256_store_pd(data, a);

	// For 4 values adding them directly is more efficient
	return data[0] + data[1] + data[2] + data[3];
}



template<>
inline typename Matrix2D_S<float>::intrinsic_type Matrix2D_S<float>::intrinsic_subtract(intrinsic_type a, intrinsic_type b) {
	return _mm256_sub_ps(a, b);
}
template<>
inline typename Matrix2D_S<double>::intrinsic_type Matrix2D_S<double>::intrinsic_subtract(intrinsic_type a, intrinsic_type b) {
	return _mm256_sub_pd(a, b);
}



template<>
inline typename Matrix2D_S<float>::intrinsic_type Matrix2D_S<float>::intrinsic_multiply(intrinsic_type a, intrinsic_type b) {
	return _mm256_mul_ps(a, b);
}
template<>
inline typename Matrix2D_S<double>::intrinsic_type Matrix2D_S<double>::intrinsic_multiply(intrinsic_type a, intrinsic_type b) {
	return _mm256_mul_pd(a, b);
}



template<>
inline typename Matrix2D_S<float>::intrinsic_type Matrix2D_S<float>::intrinsic_multiply_add(intrinsic_type a, intrinsic_type b, intrinsic_type c) {
	return _mm256_fmadd_ps(a, b, c);
}
template<>
inline typename Matrix2D_S<double>::intrinsic_type Matrix2D_S<double>::intrinsic_multiply_add(intrinsic_type a, intrinsic_type b, intrinsic_type c) {
	return _mm256_fmadd_pd(a, b, c);
}



template<>
inline typename Matrix2D_S<float>::intrinsic_type Matrix2D_S<float>::intrinsic_divide(intrinsic_type a, intrinsic_type b) {
	return _mm256_div_ps(a, b);
}
template<>
inline typename Matrix2D_S<double>::intrinsic_type Matrix2D_S<double>::intrinsic_divide(intrinsic_type a, intrinsic_type b) {
	return _mm256_div_pd(a, b);
}



template<>
inline typename Matrix2D_S<float>::intrinsic_type Matrix2D_S<float>::intrinsic_minimum(intrinsic_type a, intrinsic_type b) {
	return _mm256_min_ps(a, b);
}
template<>
inline typename Matrix2D_S<double>::intrinsic_type Matrix2D_S<double>::intrinsic_minimum(intrinsic_type a, intrinsic_type b) {
	return _mm256_min_pd(a, b);
}



template<>
inline typename Matrix2D_S<float>::intrinsic_type Matrix2D_S<float>::intrinsic_maximum(intrinsic_type a, intrinsic_type b) {
	return _mm256_max_ps(a, b);
}
template<>
inline typename Matrix2D_S<double>::intrinsic_type Matrix2D_S<double>::intrinsic_maximum(intrinsic_type a, intrinsic_type b) {
	return _mm256_max_pd(a, b);
}



template<>
inline typename Matrix2D_S<float>::intrinsic_type Matrix2D_S<float>::intrinsic_and(intrinsic_type a, intrinsic_type b) {
	return _mm256_and_ps(a, b);
}
template<>
inline typename Matrix2D_S<double>::intrinsic_type Matrix2D_S<double>::intrinsic_and(intrinsic_type a, intrinsic_type b) {
	return _mm256_and_pd(a, b);
}



template<>
inline typename Matrix2D_S<float>::intrinsic_type Matrix2D_S<float>::intrinsic_load(float* ptr) {
	return _mm256_load_ps(ptr);
}
template<>
inline typename Matrix2D_S<double>::intrinsic_type Matrix2D_S<double>::intrinsic_load(double* ptr) {
	return _mm256_load_pd(ptr);
}
template <typename T>
inline typename Matrix2D_S<T>::intrinsic_type Matrix2D_S<T>::intrinsic_load(T* ptr) {
	static_assert(std::is_same_v<T, float> || std::is_same_v<T, double>,
		"Template type must be float or double.");
	return {}; // unreachable return
}



template<>
inline typename Matrix2D_S<float>::intrinsic_type Matrix2D_S<float>::intrinsic_mask_load(mask_type* ptr) {
	return _mm256_castsi256_ps(_mm256_load_si256(reinterpret_cast<__m256i*>(ptr)));
}
template<>
inline typename Matrix2D_S<double>::intrinsic_type Matrix2D_S<double>::intrinsic_mask_load(mask_type* ptr) {
	return _mm256_castsi256_pd(_mm256_load_si256(reinterpret_cast<__m256i*>(ptr)));
}



template<>
inline void Matrix2D_S<float>::intrinsic_store(float* ptr, intrinsic_type a) {
	_mm256_store_ps(ptr, a);
}
template<>
inline void Matrix2D_S<double>::intrinsic_store(double* ptr, intrinsic_type a) {
	_mm256_store_pd(ptr, a);
}
template <typename T>
inline void Matrix2D_S<T>::intrinsic_store(T* ptr, intrinsic_type a) {
	static_assert(std::is_same_v<T, float> || std::is_same_v<T, double>,
		"Template type must be float or double.");
}



template<>
inline typename Matrix2D_S<float>::intrinsic_type Matrix2D_S<float>::intrinsic_zeroes_register() {
	return _mm256_setzero_ps();
}
template<>
inline typename Matrix2D_S<double>::intrinsic_type Matrix2D_S<double>::intrinsic_zeroes_register() {
	return _mm256_setzero_pd();
}



template<>
inline typename Matrix2D_S<float>::intrinsic_type Matrix2D_S<float>::intrinsic_init_register(float a) {
	return _mm256_set1_ps(a);
}
template<>
inline typename Matrix2D_S<double>::intrinsic_type Matrix2D_S<double>::intrinsic_init_register(double a) {
	return _mm256_set1_pd(a);
}
template <typename T>
inline typename Matrix2D_S<T>::intrinsic_type Matrix2D_S<T>::intrinsic_init_register(T a) {
	static_assert(std::is_same_v<T, float> || std::is_same_v<T, double>,
		"Template type must be float or double.");
	return {};
}



template<typename T>
T** Matrix2D_S<T>::create_H(long rows, long cols, long size, long padding) {

	// Throw an exception for a matrix size of 0.
	if (size == 0) {
		throw std::invalid_argument("Invalid matrix size, the total size (rows x columns) must be greater than 0.");
	}

	// The row pointers are created
	T** matrix = new T * [rows];

	// All the data is allocated contiguously.
	matrix[0] = static_cast<T*>(aligned_alloc(32, size * sizeof(T)));

	// Validate allocation
	if (!matrix[0]) {
		delete[] matrix;
		throw std::bad_alloc();
	}

	// Pointers are set to the contiguous block.
	for (long i = 1; i < rows; i++) {
		matrix[i] = &matrix[0][cols * i + padding * i];
	}

	return matrix;
}



template<typename T>
T** Matrix2D_S<T>::create_V(long rows, long cols, long size, long padding) {

	// Throw an exception for a matrix size of 0.
	if (size == 0) {
		throw std::invalid_argument("Invalid matrix size, the total size (rows x columns) must be greater than 0.");
	}

	// The column pointers are created
	T** matrix = new T * [cols];

	// All the data is allocated contiguously.
	matrix[0] = static_cast<T*>(aligned_alloc(32, size * sizeof(T)));

	// Validate allocation
	if (!matrix[0]) {
		delete[] matrix;
		throw std::bad_alloc();
	}

	// Pointers are set to the contiguous block.
	for (long i = 1; i < cols; i++) {
		matrix[i] = &matrix[0][rows * i + padding * i];
	}

	return matrix;
}



template<typename T>
Matrix2D_S<T>::Matrix2D_S(std::vector<std::vector<T>> data, m_orientation orientation) {

	this->rows = data.size();
	this->cols = data.at(0).size();
	this->orientation = orientation;

	// The memory is allocated to the heap.
	if (orientation == HORIZONTAL) {
		this->padding = Matrix2D_S<T>::simdStep - (this->cols % Matrix2D_S<T>::simdStep);
		this->padding = (this->padding == Matrix2D_S<T>::simdStep ? 0 : this->padding);
		this->size = (cols % Matrix2D_S<T>::simdStep == 0 ? cols : ((cols / Matrix2D_S<T>::simdStep) + 1) * Matrix2D_S<T>::simdStep) * rows;
		this->matrix = this->create_H(this->rows, this->cols, this->size, this->padding);
		this->mask = static_cast<mask_type*>(aligned_alloc(32, (this->cols + this->padding) * sizeof(T)));

		// Data is copied from a vector to an array.
		for (long i = 0; i < this->rows; i++) {

			if (this->cols != data.at(i).size()) {
				throw std::invalid_argument("Dimension inconsistencies were found in the input data.");
			}

			// Load data
			for (long j = 0; j < this->cols; j++) {
				this->matrix[i][j] = data[i][j];
			}

			// Set padding elements to 0
			for (long j = this->cols; j < this->cols + this->padding; j++) {
				this->matrix[i][j] = 0;
			}
		}

		// complete the mask
		for (long j = 0; j < this->cols + this->padding; j++) {
			this->mask[j] = -static_cast<mask_type>((j < this->cols));
		}

	} else if (orientation == VERTICAL) {
		this->padding = Matrix2D_S<T>::simdStep - (this->rows % Matrix2D_S<T>::simdStep);
		this->padding = (this->padding == Matrix2D_S<T>::simdStep ? 0 : this->padding);
		this->size = (rows % Matrix2D_S<T>::simdStep == 0 ? rows : ((rows / Matrix2D_S<T>::simdStep) + 1) * Matrix2D_S<T>::simdStep) * cols;
		this->matrix = this->create_V(this->rows, this->cols, this->size, this->padding);
		this->mask = static_cast<mask_type*>(aligned_alloc(32, (this->rows + this->padding) * sizeof(T)));

		// Validate input
		for (long i = 0; i < this->rows; i++) {
			if (this->cols != data.at(i).size()) {
				throw std::invalid_argument("Dimension inconsistencies were found in the input data.");
			}
		}

		// Data is copied from a vector to an array.
		for (long j = 0; j < this->cols; j++) {

			// Load input
			for (long i = 0; i < this->rows; i++) {
				this->matrix[j][i] = data[i][j];
			}

			// Set padding elements to 0
			for (long i = this->rows; i < this->rows + this->padding; i++) {
				this->matrix[j][i] = 0;
			}
		}

		// complete the mask
		for (long i = 0; i < this->rows + this->padding; i++) {
			this->mask[i] = -static_cast<mask_type>((i < this->rows));
		}

	} else {
		throw std::invalid_argument("No such orientation exists.");
	}
}



template<typename T>
Matrix2D_S<T>::Matrix2D_S(long rows, long cols, m_orientation orientation) {

	this->rows = rows;
	this->cols = cols;
	this->orientation = orientation;

	// The memory is allocated to the heap.
	if (orientation == HORIZONTAL) {
		this->padding = Matrix2D_S<T>::simdStep - (this->cols % Matrix2D_S<T>::simdStep);
		this->padding = (this->padding == Matrix2D_S<T>::simdStep ? 0 : this->padding);
		this->size = (cols % Matrix2D_S<T>::simdStep == 0 ? cols : ((cols / Matrix2D_S<T>::simdStep) + 1) * Matrix2D_S<T>::simdStep) * rows;
		this->matrix = this->create_H(this->rows, this->cols, this->size, this->padding);
		this->mask = static_cast<mask_type*>(aligned_alloc(32, (this->cols + this->padding) * sizeof(T)));

		// Set padding elements to 0
		for (long i = 0; i < this->rows; i++) {
			for (long j = this->cols; j < this->cols + this->padding; j++) {
				this->matrix[i][j] = 0;
			}
		}

		// complete the mask
		for (long j = 0; j < this->cols + this->padding; j++) {
			this->mask[j] = -static_cast<mask_type>((j < this->cols));
		}

	} else if (orientation == VERTICAL) {
		this->padding = Matrix2D_S<T>::simdStep - (this->rows % Matrix2D_S<T>::simdStep);
		this->padding = (this->padding == Matrix2D_S<T>::simdStep ? 0 : this->padding);
		this->size = (rows % Matrix2D_S<T>::simdStep == 0 ? rows : ((rows / Matrix2D_S<T>::simdStep) + 1) * Matrix2D_S<T>::simdStep) * cols;
		this->matrix = this->create_V(this->rows, this->cols, this->size, this->padding);
		this->mask = static_cast<mask_type*>(aligned_alloc(32, (this->rows + this->padding) * sizeof(T)));

		// Set padding elements to 0
		for (long j = 0; j < this->cols; j++) {
			for (long i = this->rows; i < this->rows + this->padding; i++) {
				this->matrix[j][i] = 0;
			}
		}

		// complete the mask
		for (long i = 0; i < this->rows + this->padding; i++) {
			this->mask[i] = -static_cast<mask_type>((i < this->rows));
		}

	} else {
		throw std::invalid_argument("No such orientation exists.");
	}
}



template<typename T>
Matrix2D_S<T>::Matrix2D_S(const Matrix2D_S& matrix) {

	this->rows = matrix.getRows();
	this->cols = matrix.getCols();
	this->orientation = matrix.getOrientation();
	T** ptr = matrix.getMatrix();

	// The memory is allocated to the heap.
	if (this->orientation == HORIZONTAL) {
		this->padding = Matrix2D_S<T>::simdStep - (this->cols % Matrix2D_S<T>::simdStep);
		this->padding = (this->padding == Matrix2D_S<T>::simdStep ? 0 : this->padding);
		this->size = (cols % Matrix2D_S<T>::simdStep == 0 ? cols : ((cols / Matrix2D_S<T>::simdStep) + 1) * Matrix2D_S<T>::simdStep) * rows;
		this->matrix = this->create_H(this->rows, this->cols, this->size, this->padding);
		this->mask = static_cast<mask_type*>(aligned_alloc(32, (this->cols + this->padding) * sizeof(T)));

		// Set padding elements to 0
		for (long i = 0; i < this->rows; i++) {
			for (long j = 0; j < this->cols + this->padding; j += Matrix2D_S<T>::simdStep) {
				Matrix2D_S<T>::intrinsic_store(&this->matrix[i][j], Matrix2D_S<T>::intrinsic_load(&ptr[i][j]));
			}
		}

		// complete the mask
		for (long j = 0; j < this->cols + this->padding; j++) {
			this->mask[j] = -static_cast<mask_type>((j < this->cols));
		}

	} else if (this->orientation == VERTICAL) {
		this->padding = Matrix2D_S<T>::simdStep - (this->rows % Matrix2D_S<T>::simdStep);
		this->padding = (this->padding == Matrix2D_S<T>::simdStep ? 0 : this->padding);
		this->size = (rows % Matrix2D_S<T>::simdStep == 0 ? rows : ((rows / Matrix2D_S<T>::simdStep) + 1) * Matrix2D_S<T>::simdStep) * cols;
		this->matrix = this->create_V(this->rows, this->cols, this->size, this->padding);
		this->mask = static_cast<mask_type*>(aligned_alloc(32, (this->rows + this->padding) * sizeof(T)));

		// Set padding elements to 0
		for (long j = 0; j < this->cols; j++) {
			for (long i = 0; i < this->rows + this->padding; i += Matrix2D_S<T>::simdStep) {
				Matrix2D_S<T>::intrinsic_store(&this->matrix[j][i], Matrix2D_S<T>::intrinsic_load(&ptr[j][i]));
			}
		}

		// complete the mask
		for (long i = 0; i < this->rows + this->padding; i++) {
			this->mask[i] = -static_cast<mask_type>((i < this->rows));
		}

	} else {
		throw std::invalid_argument("No such orientation exists.");
	}
}



template<typename T>
Matrix2D_S<T>::Matrix2D_S(long rows, long cols, m_orientation orientation, T initializer) {

	this->rows = rows;
	this->cols = cols;
	this->orientation = orientation;

	// The memory is allocated to the heap.
	if (orientation == HORIZONTAL) {
		this->padding = Matrix2D_S<T>::simdStep - (this->cols % Matrix2D_S<T>::simdStep);
		this->padding = (this->padding == Matrix2D_S<T>::simdStep ? 0 : this->padding);
		this->size = (cols % Matrix2D_S<T>::simdStep == 0 ? cols : ((cols / Matrix2D_S<T>::simdStep) + 1) * Matrix2D_S<T>::simdStep) * rows;
		this->matrix = this->create_H(this->rows, this->cols, this->size, this->padding);
		this->mask = static_cast<mask_type*>(aligned_alloc(32, (this->cols + this->padding) * sizeof(T)));

		// registers used to peform the SIMD operations
		intrinsic_type external = Matrix2D_S<T>::intrinsic_init_register(initializer);

		// Set padding elements to 0
		for (long i = 0; i < this->rows; i++) {
			for (long j = 0; j < this->cols; j += Matrix2D_S<T>::simdStep) {

				Matrix2D_S<T>::intrinsic_store(&this->matrix[i][j], external);
			}

			for (long j = this->cols; j < this->cols + this->padding; j++) {
				this->matrix[i][j] = 0;
			}
		}

		// complete the mask
		for (long j = 0; j < this->cols + this->padding; j++) {
			this->mask[j] = -static_cast<mask_type>((j < this->cols));
		}

	} else if (orientation == VERTICAL) {
		this->padding = Matrix2D_S<T>::simdStep - (this->rows % Matrix2D_S<T>::simdStep);
		this->padding = (this->padding == Matrix2D_S<T>::simdStep ? 0 : this->padding);
		this->size = (rows % Matrix2D_S<T>::simdStep == 0 ? rows : ((rows / Matrix2D_S<T>::simdStep) + 1) * Matrix2D_S<T>::simdStep) * cols;
		this->matrix = this->create_V(this->rows, this->cols, this->size, this->padding);
		this->mask = static_cast<mask_type*>(aligned_alloc(32, (this->rows + this->padding) * sizeof(T)));

		// registers used to peform the SIMD operations
		intrinsic_type external = Matrix2D_S<T>::intrinsic_init_register(initializer);

		// Set padding elements to 0
		for (long j = 0; j < this->cols; j++) {
			for (long i = 0; i < this->rows; i += Matrix2D_S<T>::simdStep) {

				Matrix2D_S<T>::intrinsic_store(&this->matrix[j][i], external);
			}

			for (long i = this->rows; i < this->rows + this->padding; i++) {
				this->matrix[j][i] = 0;
			}
		}

		// complete the mask
		for (long i = 0; i < this->rows + this->padding; i++) {
			this->mask[i] = -static_cast<mask_type>((i < this->rows));
		}

	} else {
		throw std::invalid_argument("No such orientation exists.");
	}
}



template<typename T>
Matrix2D_S<T>::Matrix2D_S(long rows, long cols, m_orientation orientation, T(*operation)()) {

	this->rows = rows;
	this->cols = cols;
	this->orientation = orientation;

	// The memory is allocated to the heap.
	if (orientation == HORIZONTAL) {
		this->padding = Matrix2D_S<T>::simdStep - (this->cols % Matrix2D_S<T>::simdStep);
		this->padding = (this->padding == Matrix2D_S<T>::simdStep ? 0 : this->padding);
		this->size = (cols % Matrix2D_S<T>::simdStep == 0 ? cols : ((cols / Matrix2D_S<T>::simdStep) + 1) * Matrix2D_S<T>::simdStep) * rows;
		this->matrix = this->create_H(this->rows, this->cols, this->size, this->padding);
		this->mask = static_cast<mask_type*>(aligned_alloc(32, (this->cols + this->padding) * sizeof(T)));

		// Set padding elements to 0
		for (long i = 0; i < this->rows; i++) {
			for (long j = 0; j < this->cols; j++) {
				this->matrix[i][j] = operation();
			}

			for (long j = this->cols; j < this->cols + this->padding; j++) {
				this->matrix[i][j] = 0;
			}
		}

		// complete the mask
		for (long j = 0; j < this->cols + this->padding; j++) {
			this->mask[j] = -static_cast<mask_type>((j < this->cols));
		}

	} else if (orientation == VERTICAL) {
		this->padding = Matrix2D_S<T>::simdStep - (this->rows % Matrix2D_S<T>::simdStep);
		this->padding = (this->padding == Matrix2D_S<T>::simdStep ? 0 : this->padding);
		this->size = (rows % Matrix2D_S<T>::simdStep == 0 ? rows : ((rows / Matrix2D_S<T>::simdStep) + 1) * Matrix2D_S<T>::simdStep) * cols;
		this->matrix = this->create_V(this->rows, this->cols, this->size, this->padding);
		this->mask = static_cast<mask_type*>(aligned_alloc(32, (this->rows + this->padding) * sizeof(T)));

		// Set padding elements to 0
		for (long j = 0; j < this->cols; j++) {
			for (long i = 0; i < this->rows; i++) {
				this->matrix[j][i] = operation();
			}

			for (long i = this->rows; i < this->rows + this->padding; i++) {
				this->matrix[j][i] = 0;
			}
		}

		// complete the mask
		for (long i = 0; i < this->rows + this->padding; i++) {
			this->mask[i] = -static_cast<mask_type>((i < this->rows));
		}

	} else {
		throw std::invalid_argument("No such orientation exists.");
	}
}



template<typename T>
Matrix2D_S<T>::~Matrix2D_S() {
	// Delete heap allocations
	aligned_free((void*)this->matrix[0]);
	aligned_free((void*)this->mask);
	delete[] this->matrix;
}



template<typename T>
Matrix2D_S<T> Matrix2D_S<T>::dot_H(Matrix2D_S<T>& matrix) {

	long numRows = matrix.getRows();
	long numCols = matrix.getCols();

	// Conditional throws exception for incompatible dimensions or orientations
	if (this->cols != numRows) {
		throw std::invalid_argument("The number of columns in matrix 1 does not equal the number of rows in matrix 2.");
	} else if (this->getOrientation() != HORIZONTAL || matrix.getOrientation() != VERTICAL) {
		throw std::invalid_argument("The orientation of matrix 1 must be horizontal, and the orientation of matrix 2 "
			"must be vertical - for simd optimization.");
	}

	T** externalMatrix = matrix.getMatrix();

	// The output matrix is created
	Matrix2D_S<T> output(this->rows, numCols, HORIZONTAL);
	T** outputMatrix = output.getMatrix();

	// registers used to peform the SIMD operations
	intrinsic_type sum;
	intrinsic_type internal;
	intrinsic_type external;

	// For each row (internal)
	for (long i = 0; i < this->rows; i++) {

		// Multiply and sum by all columns (external)
		for (long j = 0; j < numCols; j++) {

			sum = Matrix2D_S<T>::intrinsic_zeroes_register();

			// For each row (external)
			for (long k = 0; k < numRows; k += Matrix2D_S<T>::simdStep) {

				internal = Matrix2D_S<T>::intrinsic_load(&this->matrix[i][k]);
				external = Matrix2D_S<T>::intrinsic_load(&externalMatrix[j][k]);
				sum = Matrix2D_S<T>::intrinsic_multiply_add(internal, external, sum);
			}

			outputMatrix[i][j] = Matrix2D_S<T>::intrinsic_aggregate_add(sum);
		}
	}

	return output;
}



template<typename T>
Matrix2D_S<T> Matrix2D_S<T>::dot_V(Matrix2D_S<T>& matrix) {

	long numRows = matrix.getRows();
	long numCols = matrix.getCols();

	// Conditional throws exception for incompatible dimensions or orientations
	if (this->cols != numRows) {
		throw std::invalid_argument("The number of columns in matrix 1 does not equal the number of rows in matrix 2.");
	} else if (this->getOrientation() != HORIZONTAL || matrix.getOrientation() != VERTICAL) {
		throw std::invalid_argument("The orientation of matrix 1 must be horizontal, and the orientation of matrix 2 "
			"must be vertical - for simd optimization.");
	}

	T** externalMatrix = matrix.getMatrix();

	// The output matrix is created
	Matrix2D_S<T> output(this->rows, numCols, VERTICAL);
	T** outputMatrix = output.getMatrix();

	// registers used to peform the SIMD operations
	intrinsic_type sum;
	intrinsic_type internal;
	intrinsic_type external;

	// For each row (internal)
	for (long i = 0; i < this->rows; i++) {

		// Multiply and sum by all columns (external)
		for (long j = 0; j < numCols; j++) {

			sum = Matrix2D_S<T>::intrinsic_zeroes_register();

			// For each row (external)
			for (long k = 0; k < numRows; k += Matrix2D_S<T>::simdStep) {

				internal = Matrix2D_S<T>::intrinsic_load(&this->matrix[i][k]);
				external = Matrix2D_S<T>::intrinsic_load(&externalMatrix[j][k]);
				sum = Matrix2D_S<T>::intrinsic_multiply_add(internal, external, sum);

			}

			outputMatrix[j][i] = Matrix2D_S<T>::intrinsic_aggregate_add(sum);
		}
	}

	return output;
}



template<typename T>
Matrix2D_S<T> Matrix2D_S<T>::multiply(T scalar) {

	// The output matrix is created
	Matrix2D_S<T> output(this->rows, this->cols, this->orientation);
	T** outputMatrix = output.getMatrix();

	// registers used to peform the SIMD operations
	intrinsic_type internal;
	intrinsic_type external = Matrix2D_S<T>::intrinsic_init_register(scalar);

	// Multiply each element by the scalar
	for (long i = 0; i < this->size; i += Matrix2D_S<T>::simdStep) {
		internal = Matrix2D_S<T>::intrinsic_load(&this->matrix[0][i]);
		Matrix2D_S<T>::intrinsic_store(&outputMatrix[0][i], Matrix2D_S<T>::intrinsic_multiply(internal, external));
	}

	return output;
}



template<typename T>
Matrix2D_S<T> Matrix2D_S<T>::multiply(Matrix2D_S<T>& matrix) {

	// Validate input matrix
	if (matrix.getCols() != this->cols || matrix.getRows() != this->rows || matrix.getOrientation() != this->orientation) {
		throw std::invalid_argument("The input matrix is incompatible either in it's dimensions or orientation.");
	}

	// Input matrix
	T** inputMatrix = matrix.getMatrix();

	// The output matrix is created
	Matrix2D_S<T> output(this->rows, this->cols, this->orientation);
	T** outputMatrix = output.getMatrix();

	// registers used to peform the SIMD operations
	intrinsic_type internal;
	intrinsic_type external;

	// Multiply each element by the scalar
	for (long i = 0; i < this->size; i += Matrix2D_S<T>::simdStep) {
		internal = Matrix2D_S<T>::intrinsic_load(&this->matrix[0][i]);
		external = Matrix2D_S<T>::intrinsic_load(&inputMatrix[0][i]);
		Matrix2D_S<T>::intrinsic_store(&outputMatrix[0][i], Matrix2D_S<T>::intrinsic_multiply(internal, external));
	}

	return output;
}



template<typename T>
Matrix2D_S<T> Matrix2D_S<T>::divide(T scalar) {

	// The output matrix is created
	Matrix2D_S<T> output(this->rows, this->cols, this->orientation);
	T** outputMatrix = output.getMatrix();

	// calculate the reciprocal
	T reciprocal = 1 / scalar;

	// registers used to peform the SIMD operations
	intrinsic_type internal;
	intrinsic_type external = Matrix2D_S<T>::intrinsic_init_register(reciprocal);

	// Multiply each element by the scalar
	for (long i = 0; i < this->size; i += Matrix2D_S<T>::simdStep) {
		internal = Matrix2D_S<T>::intrinsic_load(&this->matrix[0][i]);
		Matrix2D_S<T>::intrinsic_store(&outputMatrix[0][i], Matrix2D_S<T>::intrinsic_multiply(internal, external));
	}

	return output;
}



template<typename T>
Matrix2D_S<T> Matrix2D_S<T>::divide(Matrix2D_S<T>& matrix) {

	// Validate input matrix
	if (matrix.getCols() != this->cols || matrix.getRows() != this->rows || matrix.getOrientation() != this->orientation) {
		throw std::invalid_argument("The input matrix is incompatible either in it's dimensions or orientation.");
	}

	// Input matrix
	T** inputMatrix = matrix.getMatrix();

	// The output matrix is created
	Matrix2D_S<T> output(this->rows, this->cols, this->orientation);
	T** outputMatrix = output.getMatrix();

	// registers used to peform the SIMD operations
	intrinsic_type internal;
	intrinsic_type external;

	// Divide each element by the corresponding element in the input
	for (long i = 0; i < this->size; i += Matrix2D_S<T>::simdStep) {
		internal = Matrix2D_S<T>::intrinsic_load(&this->matrix[0][i]);
		external = Matrix2D_S<T>::intrinsic_load(&inputMatrix[0][i]);
		Matrix2D_S<T>::intrinsic_store(&outputMatrix[0][i], Matrix2D_S<T>::intrinsic_divide(internal, external));
	}

	return output;
}



template<typename T>
Matrix2D_S<T> Matrix2D_S<T>::add(T scalar) {

	// The output matrix is created
	Matrix2D_S<T> output(this->rows, this->cols, this->orientation);
	T** outputMatrix = output.getMatrix();

	// registers used to peform the SIMD operations
	intrinsic_type internal;
	intrinsic_type external = Matrix2D_S<T>::intrinsic_init_register(scalar);

	// The scalar is added to each element in the matrix
	if (this->orientation == HORIZONTAL) {
		for (long i = 0; i < this->rows; i++) {
			for (long j = 0; j < this->cols; j += Matrix2D_S<T>::simdStep) {
				internal = Matrix2D_S<T>::intrinsic_load(&this->matrix[i][j]);
				Matrix2D_S<T>::intrinsic_store(&outputMatrix[i][j], Matrix2D_S<T>::intrinsic_add(internal, external));
			}
		}
	} else {
		for (long j = 0; j < this->cols; j++) {
			for (long i = 0; i < this->rows; i += Matrix2D_S<T>::simdStep) {
				internal = Matrix2D_S<T>::intrinsic_load(&this->matrix[j][i]);
				Matrix2D_S<T>::intrinsic_store(&outputMatrix[j][i], Matrix2D_S<T>::intrinsic_add(internal, external));
			}
		}
	}

	return output;
}



template<typename T>
Matrix2D_S<T> Matrix2D_S<T>::add(Matrix2D_S<T>& matrix) {

	// Validate input matrix
	if (matrix.getCols() != this->cols || matrix.getRows() != this->rows || matrix.getOrientation() != this->orientation) {
		throw std::invalid_argument("The input matrix is incompatible either in it's dimensions or orientation.");
	}

	// Input matrix
	T** inputMatrix = matrix.getMatrix();

	// The output matrix is created
	Matrix2D_S<T> output(this->rows, this->cols, this->orientation);
	T** outputMatrix = output.getMatrix();

	// registers used to peform the SIMD operations
	intrinsic_type internal;
	intrinsic_type external;

	// Add each element to the corresponding element in the input matrix
	for (long i = 0; i < this->size; i += Matrix2D_S<T>::simdStep) {
		internal = Matrix2D_S<T>::intrinsic_load(&this->matrix[0][i]);
		external = Matrix2D_S<T>::intrinsic_load(&inputMatrix[0][i]);
		Matrix2D_S<T>::intrinsic_store(&outputMatrix[0][i], Matrix2D_S<T>::intrinsic_add(internal, external));
	}

	return output;
}



template<typename T>
Matrix2D_S<T> Matrix2D_S<T>::subtract(T scalar) {

	// The output matrix is created
	Matrix2D_S<T> output(this->rows, this->cols, this->orientation);
	T** outputMatrix = output.getMatrix();

	// registers used to peform the SIMD operations
	intrinsic_type internal;
	intrinsic_type external = Matrix2D_S<T>::intrinsic_init_register(scalar);

	// The scalar is subtracted from every element in the matrix
	if (this->orientation == HORIZONTAL) {
		for (long i = 0; i < this->rows; i++) {
			for (long j = 0; j < this->cols; j += Matrix2D_S<T>::simdStep) {
				internal = Matrix2D_S<T>::intrinsic_load(&this->matrix[i][j]);
				Matrix2D_S<T>::intrinsic_store(&outputMatrix[i][j], Matrix2D_S<T>::intrinsic_subtract(internal, external));
			}
		}
	} else {
		for (long j = 0; j < this->cols; j++) {
			for (long i = 0; i < this->rows; i += Matrix2D_S<T>::simdStep) {
				internal = Matrix2D_S<T>::intrinsic_load(&this->matrix[j][i]);
				Matrix2D_S<T>::intrinsic_store(&outputMatrix[j][i], Matrix2D_S<T>::intrinsic_subtract(internal, external));
			}
		}
	}

	return output;
}



template<typename T>
Matrix2D_S<T> Matrix2D_S<T>::subtract(Matrix2D_S<T>& matrix) {

	// Validate input matrix
	if (matrix.getCols() != this->cols || matrix.getRows() != this->rows || matrix.getOrientation() != this->orientation) {
		throw std::invalid_argument("The input matrix is incompatible either in it's dimensions or orientation.");
	}

	// Input matrix
	T** inputMatrix = matrix.getMatrix();

	// The output matrix is created
	Matrix2D_S<T> output(this->rows, this->cols, this->orientation);
	T** outputMatrix = output.getMatrix();

	// registers used to peform the SIMD operations
	intrinsic_type internal;
	intrinsic_type external;

	// Subtract the corresponding element in the input matrix from a given element
	for (long i = 0; i < this->size; i += Matrix2D_S<T>::simdStep) {
		internal = Matrix2D_S<T>::intrinsic_load(&this->matrix[0][i]);
		external = Matrix2D_S<T>::intrinsic_load(&inputMatrix[0][i]);
		Matrix2D_S<T>::intrinsic_store(&outputMatrix[0][i], Matrix2D_S<T>::intrinsic_subtract(internal, external));
	}

	return output;
}



template<typename T>
Matrix2D_S<T> Matrix2D_S<T>::min(T scalar) {

	// The output matrix is created
	Matrix2D_S<T> output(this->rows, this->cols, this->orientation);
	T** outputMatrix = output.getMatrix();

	// registers used to peform the SIMD operations
	intrinsic_type internal;
	intrinsic_type external = Matrix2D_S<T>::intrinsic_init_register(scalar);
	intrinsic_type mask;

	// The minimum is found relative to the scalar
	if (this->orientation == HORIZONTAL) {
		for (long i = 0; i < this->rows; i++) {
			for (long j = 0; j < this->cols; j += Matrix2D_S<T>::simdStep) {
				mask = Matrix2D_S<T>::intrinsic_mask_load(&this->mask[j]);
				internal = Matrix2D_S<T>::intrinsic_load(&this->matrix[i][j]);
				Matrix2D_S<T>::intrinsic_store(&outputMatrix[i][j], Matrix2D_S<T>::intrinsic_and(Matrix2D_S<T>::intrinsic_minimum(internal, external), mask));
			}
		}
	} else {
		for (long j = 0; j < this->cols; j++) {
			for (long i = 0; i < this->rows; i += Matrix2D_S<T>::simdStep) {
				mask = Matrix2D_S<T>::intrinsic_mask_load(&this->mask[i]);
				internal = Matrix2D_S<T>::intrinsic_load(&this->matrix[j][i]);
				Matrix2D_S<T>::intrinsic_store(&outputMatrix[j][i], Matrix2D_S<T>::intrinsic_and(Matrix2D_S<T>::intrinsic_minimum(internal, external), mask));
			}
		}
	}

	return output;
}



template<typename T>
Matrix2D_S<T> Matrix2D_S<T>::min(Matrix2D_S<T>& matrix) {

	// Validate input matrix
	if (matrix.getCols() != this->cols || matrix.getRows() != this->rows || matrix.getOrientation() != this->orientation) {
		throw std::invalid_argument("The input matrix is incompatible either in it's dimensions or orientation.");
	}

	// Input matrix
	T** inputMatrix = matrix.getMatrix();

	// The output matrix is created
	Matrix2D_S<T> output(this->rows, this->cols, this->orientation);
	T** outputMatrix = output.getMatrix();

	// registers used to peform the SIMD operations
	intrinsic_type internal;
	intrinsic_type external;

	// Find the elementwise minimum of two matrices 
	for (long i = 0; i < this->size; i += Matrix2D_S<T>::simdStep) {
		internal = Matrix2D_S<T>::intrinsic_load(&this->matrix[0][i]);
		external = Matrix2D_S<T>::intrinsic_load(&inputMatrix[0][i]);
		Matrix2D_S<T>::intrinsic_store(&outputMatrix[0][i], Matrix2D_S<T>::intrinsic_minimum(internal, external));
	}

	return output;
}



template<typename T>
Matrix2D_S<T> Matrix2D_S<T>::max(T scalar) {

	// The output matrix is created
	Matrix2D_S<T> output(this->rows, this->cols, this->orientation);
	T** outputMatrix = output.getMatrix();

	// registers used to peform the SIMD operations
	intrinsic_type internal;
	intrinsic_type external = Matrix2D_S<T>::intrinsic_init_register(scalar);
	intrinsic_type mask;

	// The maximum value is found relative to the scalar
	if (this->orientation == HORIZONTAL) {
		for (long i = 0; i < this->rows; i++) {
			for (long j = 0; j < this->cols; j += Matrix2D_S<T>::simdStep) {
				mask = Matrix2D_S<T>::intrinsic_mask_load(&this->mask[j]);
				internal = Matrix2D_S<T>::intrinsic_load(&this->matrix[i][j]);
				Matrix2D_S<T>::intrinsic_store(&outputMatrix[i][j], Matrix2D_S<T>::intrinsic_and(Matrix2D_S<T>::intrinsic_maximum(internal, external), mask));
			}
		}
	} else {
		for (long j = 0; j < this->cols; j++) {
			for (long i = 0; i < this->rows; i += Matrix2D_S<T>::simdStep) {
				mask = Matrix2D_S<T>::intrinsic_mask_load(&this->mask[i]);
				internal = Matrix2D_S<T>::intrinsic_load(&this->matrix[j][i]);
				Matrix2D_S<T>::intrinsic_store(&outputMatrix[j][i], Matrix2D_S<T>::intrinsic_and(Matrix2D_S<T>::intrinsic_maximum(internal, external), mask));
			}
		}
	}

	return output;
}



template<typename T>
Matrix2D_S<T> Matrix2D_S<T>::max(Matrix2D_S<T>& matrix) {

	// Validate input matrix
	if (matrix.getCols() != this->cols || matrix.getRows() != this->rows || matrix.getOrientation() != this->orientation) {
		throw std::invalid_argument("The input matrix is incompatible either in it's dimensions or orientation.");
	}

	// Input matrix
	T** inputMatrix = matrix.getMatrix();

	// The output matrix is created
	Matrix2D_S<T> output(this->rows, this->cols, this->orientation);
	T** outputMatrix = output.getMatrix();

	// registers used to peform the SIMD operations
	intrinsic_type internal;
	intrinsic_type external;

	// Find the elementwise maximum of two matrices 
	for (long i = 0; i < this->size; i += Matrix2D_S<T>::simdStep) {
		internal = Matrix2D_S<T>::intrinsic_load(&this->matrix[0][i]);
		external = Matrix2D_S<T>::intrinsic_load(&inputMatrix[0][i]);
		Matrix2D_S<T>::intrinsic_store(&outputMatrix[0][i], Matrix2D_S<T>::intrinsic_maximum(internal, external));
	}

	return output;
}



template<typename T>
Matrix2D_S<T> Matrix2D_S<T>::apply(T(*operation)(T)) {

	// The output matrix is created
	Matrix2D_S<T> output(this->rows, this->cols, this->orientation);
	T** outputMatrix = output.getMatrix();

	// The operation is performed
	if (this->orientation == HORIZONTAL) {
		for (long i = 0; i < this->rows; i++) {
			for (long j = 0; j < this->cols; j++) {
				outputMatrix[i][j] = operation(this->matrix[i][j]);
			}
		}
	} else {
		for (long j = 0; j < this->cols; j++) {
			for (long i = 0; i < this->rows; i++) {
				outputMatrix[j][i] = operation(this->matrix[j][i]);
			}
		}
	}

	return output;
}



template<typename T>
Matrix2D_S<T> Matrix2D_S<T>::transpose() {

	// The output matrix is created
	Matrix2D_S<T> output(this->cols, this->rows, this->orientation);
	T** outputMatrix = output.getMatrix();

	// The matrix is transposed
	if (this->orientation == HORIZONTAL) {
		for (long i = 0; i < this->rows; i++) {
			for (long j = 0; j < this->cols; j++) {
				outputMatrix[j][i] = this->matrix[i][j];
			}
		}
	} else {
		for (long j = 0; j < this->cols; j++) {
			for (long i = 0; i < this->rows; i++) {
				outputMatrix[i][j] = this->matrix[j][i];
			}
		}
	}

	return output;
}



template<typename T>
Matrix2D_S<T> Matrix2D_S<T>::realign() {

	if (this->orientation == HORIZONTAL) {

		Matrix2D_S<T> output(this->rows, this->cols, VERTICAL);
		T** outputMatrix = output.getMatrix();

		// re-align the memory
		for (long i = 0; i < this->rows; i++) {
			for (long j = 0; j < this->cols; j++) {
				outputMatrix[j][i] = this->matrix[i][j];
			}
		}

		return output;

	} else {

		Matrix2D_S<T> output(this->rows, this->cols, HORIZONTAL);
		T** outputMatrix = output.getMatrix();

		// re-align the memory
		for (long i = 0; i < this->rows; i++) {
			for (long j = 0; j < this->cols; j++) {
				outputMatrix[i][j] = this->matrix[j][i];
			}
		}

		return output;
	}
}



template<typename T>
long Matrix2D_S<T>::getRows() {
	return this->rows;
}



template<typename T>
long Matrix2D_S<T>::getCols() {
	return this->cols;
}



template<typename T>
long Matrix2D_S<T>::getPadding() {
	return this->padding;
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

	if (orientation == HORIZONTAL) {

		// Iterate over matrix and convert to string data.
		for (long i = 0; i < this->rows; i++) {
			for (long j = 0; j < this->cols; j++) {
				matrix += std::to_string(this->matrix[i][j]) + ", ";
			}
			matrix += "\n";
		}

	} else {

		// Iterate over matrix and convert to string data.
		for (long i = 0; i < this->rows; i++) {
			for (long j = 0; j < this->cols; j++) {
				matrix += std::to_string(this->matrix[j][i]) + ", ";
			}
			matrix += "\n";
		}

	}

	return matrix;
}