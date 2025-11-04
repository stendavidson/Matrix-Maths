#pragma once

#include <vector>
#include <stdexcept>
#include <string>
#include <sstream>
#include <format>
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



// The memory allocation function(s) and compiler specific features are defined
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

	/**
	 * The RESTRICT keyword lets the compiler know that the same memory is not being
	 * accessed via two different pointers. This is the MSVC version.
	 */
	#define RESTRICT __restrict


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

	/**
	 * The RESTRICT keyword lets the compiler know that the same memory is not being
	 * accessed via two different pointers. This is the GCC/Clang version.
	 */
	#define RESTRICT __restrict__

#endif


// Define alignments
enum class M_Orientation: int { HORIZONTAL, VERTICAL };


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
	Matrix2D_S(const std::vector<std::vector<T>> data, const M_Orientation orientation);

	/**
	 * Constructor - creates a matrix of a given size.
	 *
	 * @param rows: The number of rows in the matrix being instantiated.
	 *
	 * @param cols: The number of columsn int the matrix being instantiated.
	 *
	 * @param orientation: The memory orientation of the matrix being created.
	 */
	Matrix2D_S(const size_t rows, const size_t cols, const M_Orientation orientation);

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
	Matrix2D_S(const size_t rows, const size_t cols, const M_Orientation orientation, const T initializer);

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
	Matrix2D_S(const size_t rows, const size_t cols, const M_Orientation orientation, T(* const operation)());

	/**
	 * Destructor
	 */
	~Matrix2D_S();

	/**
	 * This method finds the dot product of two matrices - output has M_Orientation::HORIZONTALly
	 * aligned memory.
	 *
	 * @param matrix: The matrix with which the dot product is being performed.
	 */
	Matrix2D_S<T> dot_H(const Matrix2D_S<T>& matrix) const;

	/**
	* This method finds the dot product of two matrices - output has M_Orientation::VERTICALly
	* aligned memory.
	*
	* @param matrix: The matrix with which the dot product is being performed.
	*/
	Matrix2D_S<T> dot_V(const Matrix2D_S<T>& matrix) const;

	/**
	* This method finds the multiplication of a salar and a matrice
	*
	* @param scalar: The scalar with which to multiply this matrix.
	*/
	Matrix2D_S<T> multiply(const T scalar) const;

	/**
	* This method finds the elementwise multiplication of two matrices.
	*
	* @param scalar: The matrix with which to multiply this matrix.
	*/
	Matrix2D_S<T> multiply(const Matrix2D_S<T>& matrix) const;

	/**
	* This method finds the multiplication of a salar and a matrice
	*
	* @param scalar: The scalar with which to divide this matrix.
	*/
	Matrix2D_S<T> divide(const T scalar) const;

	/**
	* This method finds the elementwise division of two matrices.
	*
	* @param scalar: The matrix with which to divide this matrix.
	*/
	Matrix2D_S<T> divide(const Matrix2D_S<T>& matrix) const;

	/**
	* This method finds the addition of a scalar to a matrice
	*
	* @param scalar: The scalar with which to add this matrix.
	*/
	Matrix2D_S<T> add(const T scalar) const;

	/**
	* This method finds the addition of two matrices
	*
	* @param matrix: The matrix with which to add this matrix.
	*/
	Matrix2D_S<T> add(const Matrix2D_S<T>& matrix) const;

	/**
	* This method finds the subtaction of a scalar to a matrice
	*
	* @param scalar: The scalar with which to subract this matrix.
	*/
	Matrix2D_S<T> subtract(const T scalar) const;

	/**
	* This method finds the subtraction of a matrice.
	*
	* @param matrix: The matrix with which to subtract this matrix.
	*/
	Matrix2D_S<T> subtract(const Matrix2D_S<T>& matrix) const;

	/**
	* This method finds the elementwise max of a scalar.
	*
	* @param scalar: The scalar against which the maximum is found.
	*/
	Matrix2D_S<T> max(const T scalar) const;

	/**
	* This method finds the maximum values of two matrices.
	*
	* @param matrix: The matrix against which to find the maximum
	*/
	Matrix2D_S<T> max(const Matrix2D_S<T>& matrix) const;

	/**
	* This method finds the elementwise min of a scalar
	*
	* @param scalar: The scalar against which the minimum is found.
	*/
	Matrix2D_S<T> min(const T scalar) const;

	/**
	* This method finds the minimum values of two matrices.
	*
	* @param matrix: The matrix against which to find the minimum
	*/
	Matrix2D_S<T> min(const Matrix2D_S<T>& matrix) const;

	/**
	* This method applies a single operation to the entire matrix elementwise.
	*
	* @param operation: The operation which is performed on the entire array.
	*/
	Matrix2D_S<T> apply(T(* const operation)(T)) const;

	/**
	* This method returns the Matrix2D_S transpose of the given matrice.
	*/
	Matrix2D_S<T> transpose() const;

	/**
	* This method re-allocates the memory so that it is aligned perpendicularly
	* to how it was previously.
	*/
	Matrix2D_S<T> realign() const;
	
	/**
	 * This access operator allows the underlying data to be more easily and safely accessed
	 */
	T& operator()(const size_t row, const size_t col);

	/**
	 * This access operator allows the underlying data to be more easily and safely accessed
	 */
	T& operator()(const size_t row, const size_t col) const;

	/**
	 * This getter method returns the the number of rows in the matrix
	 */
	size_t getRows() const noexcept;

	/**
	 * This getter method returns the the number of rows in the matrix
	 */
	size_t getCols() const noexcept;

	/**
	 * This getter method returns the amount of padding per row or column (depending on memory orientation)
	 */
	size_t getPadding() const noexcept;

	/**
	 * This getter method returns the memory orientation of the array.
	 */
	M_Orientation getOrientation() const noexcept;

	/**
	 * This method converts the Matrix to a string.
	 */
	std::string toString() const noexcept;

private:

	/**
	 * The number of rows.
	 */
	size_t rows;

	/**
	 * The number of columns.
	 */
	size_t cols;

	/**
	 * The padded dimension
	 */
	size_t paddedDim;

	/**
	 * The total size of the array including padding
	 */
	size_t size;

	/**
	 * The column-wise or row-wise padding
	 */
	size_t padding;

	/**
	 * The orientation of memory as assigned for the Matrix.
	 */
	M_Orientation orientation;

	/**
	 * This nested array holds the matrix data.
	 */
	T* RESTRICT matrix;

	/**
	 * This array holds the row or column mask
	 */
	mask_type* mask;

	/**
	 * This method performs creates a new memory aligned contiguously allocated array.
	 *
	 * @param size: The total number of elements in the array.
	 */
	T* create(const size_t size);

	/**
	 * This static value holds the SIMD step size.
	 */
	static size_t simdStep;

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
	 * This static value holds the maximum value possible for type T
	 */
	static T maximum;

	/**
	 * This static value holds the minimum value possible for type T
	 */
	static T minimum;
};



template<typename T>
T Matrix2D_S<T>::minimum = std::numeric_limits<T>::min();



template<typename T>
T Matrix2D_S<T>::maximum = std::numeric_limits<T>::max();



template<>
size_t Matrix2D_S<float>::simdStep = 8;
template<>
size_t Matrix2D_S<double>::simdStep = 4;



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

	// Array aggregated M_Orientation::HORIZONTALly - for 8 values this is sees a performance improvement
	__m128 half1 = _mm256_castps256_ps128(a);
	__m128 half2 = _mm256_extractf128_ps(a, 1);
	__m128 sum128 = _mm_add_ps(half1, half2);

	// Now M_Orientation::HORIZONTALly add the 4 floats in sum128.
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
T* Matrix2D_S<T>::create(const size_t size) {

	// Throw an exception for a matrix size of 0.
	if (size == 0) {
		throw std::invalid_argument("Invalid matrix size, the total size (rows x columns) must be greater than 0.");
	}

	// The matrix memory is allocated
	T* matrix = static_cast<T*>(aligned_alloc(32, size * sizeof(T)));

	// Validate allocation
	if (!matrix) {
		delete[] matrix;
		throw std::bad_alloc();
	}

	return matrix;
}



template<typename T>
Matrix2D_S<T>::Matrix2D_S(const std::vector<std::vector<T>> data, const M_Orientation orientation) {

	this->rows = data.size();
	this->cols = data.at(0).size();
	this->orientation = orientation;

	// The memory is allocated to the heap.
	if (orientation == M_Orientation::HORIZONTAL) {
		this->padding = Matrix2D_S<T>::simdStep - (this->cols % Matrix2D_S<T>::simdStep);
		this->padding = (this->padding == Matrix2D_S<T>::simdStep ? 0 : this->padding);
		this->paddedDim = this->cols + this->padding;
		this->size = this->paddedDim * this->rows;
		this->matrix = this->create(this->size);
		this->mask = static_cast<mask_type*>(aligned_alloc(32, this->paddedDim * sizeof(T)));


		// Offset used to accumulate instead of multiplying
		size_t dimOffset = 0;

		// Data is copied from a vector to an array.
		for (size_t i = 0; i < this->rows; i++) {

			if (this->cols != data.at(i).size()) {
				throw std::invalid_argument(std::format("Row {} has an inconsistent number of values.", i + 1));
			}

			// Load data
			for (size_t j = 0; j < this->cols; j++) {
				this->matrix[dimOffset + j] = data[i][j];
			}

			// Set padding elements to 0
			for (size_t j = this->cols; j < this->paddedDim; j++) {
				this->matrix[dimOffset + j] = 0;
			}

			// Accumulate each row
			dimOffset += this->paddedDim;
		}

		// complete the mask
		for (size_t j = 0; j < this->paddedDim; j++) {
			this->mask[j] = -static_cast<mask_type>((j < this->cols));
		}

	} else if (orientation == M_Orientation::VERTICAL) {
		this->padding = Matrix2D_S<T>::simdStep - (this->rows % Matrix2D_S<T>::simdStep);
		this->padding = (this->padding == Matrix2D_S<T>::simdStep ? 0 : this->padding);
		this->paddedDim = this->rows + this->padding;
		this->size = this->paddedDim * this->cols;
		this->matrix = this->create(this->size);
		this->mask = static_cast<mask_type*>(aligned_alloc(32, this->paddedDim * sizeof(T)));

		// Validate input
		for (size_t i = 0; i < this->rows; i++) {
			if(this->cols != data.at(i).size()) {
				throw std::invalid_argument(std::format("Row {} has an inconsistent number of values.", i + 1));
			}
		}

		// Offset used to accumulate instead of multiplying
		size_t dimOffset = 0;

		// Data is copied from a vector to an array.
		for (size_t j = 0; j < this->cols; j++) {

			// Load input
			for (size_t i = 0; i < this->rows; i++) {
				this->matrix[dimOffset + i] = data[i][j];
			}

			// Set padding elements to 0
			for (size_t i = this->rows; i < this->paddedDim; i++) {
				this->matrix[dimOffset + i] = 0;
			}

			// Accumulate each column
			dimOffset += this->paddedDim;
		}

		// complete the mask
		for (size_t i = 0; i < this->paddedDim; i++) {
			this->mask[i] = -static_cast<mask_type>((i < this->rows));
		}

	}
}



template<typename T>
Matrix2D_S<T>::Matrix2D_S(const size_t rows, const size_t cols, M_Orientation orientation) {

	this->rows = rows;
	this->cols = cols;
	this->orientation = orientation;

	// The memory is allocated to the heap.
	if (orientation == M_Orientation::HORIZONTAL) {
		this->padding = Matrix2D_S<T>::simdStep - (this->cols % Matrix2D_S<T>::simdStep);
		this->padding = (this->padding == Matrix2D_S<T>::simdStep ? 0 : this->padding);
		this->paddedDim = this->cols + this->padding;
		this->size = this->paddedDim * this->rows;
		this->matrix = this->create(this->size);
		this->mask = static_cast<mask_type*>(aligned_alloc(32, this->paddedDim * sizeof(T)));

		// Offset used to accumulate instead of multiplying
		size_t dimOffset = 0;

		// Set padding elements to 0
		for (size_t i = 0; i < this->rows; i++) {
			
			for (size_t j = this->cols; j < this->paddedDim; j++) {
				this->matrix[dimOffset + j] = 0;
			}

			// Accumulate each row
			dimOffset += this->paddedDim;
		}

		// complete the mask
		for (size_t j = 0; j < this->paddedDim; j++) {
			this->mask[j] = -static_cast<mask_type>((j < this->cols));
		}

	} else if (orientation == M_Orientation::VERTICAL) {
		this->padding = Matrix2D_S<T>::simdStep - (this->rows % Matrix2D_S<T>::simdStep);
		this->padding = (this->padding == Matrix2D_S<T>::simdStep ? 0 : this->padding);
		this->paddedDim = this->rows + this->padding;
		this->size = this->paddedDim * this->cols;
		this->matrix = this->create(this->size);
		this->mask = static_cast<mask_type*>(aligned_alloc(32, this->paddedDim * sizeof(T)));

		// Offset used to accumulate instead of multiplying
		size_t dimOffset = 0;

		// Set padding elements to 0
		for (size_t j = 0; j < this->cols; j++) {
			
			for (size_t i = this->rows; i < this->paddedDim; i++) {
				this->matrix[dimOffset + i] = 0;
			}

			// Accumulate each column
			dimOffset += this->paddedDim;
		}

		// complete the mask
		for (size_t i = 0; i < this->paddedDim; i++) {
			this->mask[i] = -static_cast<mask_type>((i < this->rows));
		}

	}
}



template<typename T>
Matrix2D_S<T>::Matrix2D_S(const Matrix2D_S& matrix) {

	this->rows = matrix.getRows();
	this->cols = matrix.getCols();
	this->orientation = matrix.getOrientation();

	// The memory is allocated to the heap.
	if (this->orientation == M_Orientation::HORIZONTAL) {
		this->padding = Matrix2D_S<T>::simdStep - (this->cols % Matrix2D_S<T>::simdStep);
		this->padding = (this->padding == Matrix2D_S<T>::simdStep ? 0 : this->padding);
		this->paddedDim = this->cols + this->padding;
		this->size = this->paddedDim * this->rows;
		this->matrix = this->create(this->size);
		this->mask = static_cast<mask_type*>(aligned_alloc(32, this->paddedDim * sizeof(T)));

		// Offset used to accumulate instead of multiplying
		size_t dimOffset = 0;

		// Set padding elements to 0
		for (size_t i = 0; i < this->rows; i++) {

			for (size_t j = 0; j < this->paddedDim; j += Matrix2D_S<T>::simdStep) {
				Matrix2D_S<T>::intrinsic_store(&this->matrix[dimOffset + j], Matrix2D_S<T>::intrinsic_load(&matrix.matrix[dimOffset + j]));
			}

			// Accumulate each row
			dimOffset += this->paddedDim;
		}

		// complete the mask
		for (size_t j = 0; j < this->paddedDim; j++) {
			this->mask[j] = -static_cast<mask_type>((j < this->cols));
		}

	} else if (this->orientation == M_Orientation::VERTICAL) {
		this->padding = Matrix2D_S<T>::simdStep - (this->rows % Matrix2D_S<T>::simdStep);
		this->padding = (this->padding == Matrix2D_S<T>::simdStep ? 0 : this->padding);
		this->paddedDim = this->rows + this->padding;
		this->size = this->paddedDim * this->cols;
		this->matrix = this->create(this->size);
		this->mask = static_cast<mask_type*>(aligned_alloc(32, this->paddedDim * sizeof(T)));

		// Offset used to accumulate instead of multiplying
		size_t dimOffset = 0;

		// Set padding elements to 0
		for (size_t j = 0; j < this->cols; j++) {

			for (size_t i = 0; i < this->paddedDim; i += Matrix2D_S<T>::simdStep) {
				Matrix2D_S<T>::intrinsic_store(&this->matrix[dimOffset + i], Matrix2D_S<T>::intrinsic_load(&matrix.matrix[dimOffset + i]));
			}

			// Accumulate each column
			dimOffset += this->paddedDim;
		}

		// complete the mask
		for (size_t i = 0; i < this->paddedDim; i++) {
			this->mask[i] = -static_cast<mask_type>((i < this->rows));
		}

	}
}



template<typename T>
Matrix2D_S<T>::Matrix2D_S(const size_t rows, const size_t cols, const M_Orientation orientation, const T initializer) {

	this->rows = rows;
	this->cols = cols;
	this->orientation = orientation;

	// The memory is allocated to the heap.
	if (orientation == M_Orientation::HORIZONTAL) {
		this->padding = Matrix2D_S<T>::simdStep - (this->cols % Matrix2D_S<T>::simdStep);
		this->padding = (this->padding == Matrix2D_S<T>::simdStep ? 0 : this->padding);
		this->paddedDim = this->cols + this->padding;
		this->size = this->paddedDim * this->rows;
		this->matrix = this->create(this->size);
		this->mask = static_cast<mask_type*>(aligned_alloc(32, this->paddedDim * sizeof(T)));

		// Offset used to accumulate instead of multiplying
		size_t dimOffset = 0;

		// registers used to peform the SIMD operations
		intrinsic_type external = Matrix2D_S<T>::intrinsic_init_register(initializer);

		// Set padding elements to 0
		for (size_t i = 0; i < this->rows; i++) {
			
			for (size_t j = 0; j < this->cols; j += Matrix2D_S<T>::simdStep) {
				Matrix2D_S<T>::intrinsic_store(&this->matrix[dimOffset + j], external);
			}

			for (size_t j = this->cols; j < this->paddedDim; j++) {
				this->matrix[dimOffset + j] = 0;
			}

			// Accumulate each row
			dimOffset += this->paddedDim;
		}

		// complete the mask
		for (size_t j = 0; j < this->paddedDim; j++) {
			this->mask[j] = -static_cast<mask_type>((j < this->cols));
		}

	} else if (orientation == M_Orientation::VERTICAL) {
		this->padding = Matrix2D_S<T>::simdStep - (this->rows % Matrix2D_S<T>::simdStep);
		this->padding = (this->padding == Matrix2D_S<T>::simdStep ? 0 : this->padding);
		this->paddedDim = this->rows + this->padding;
		this->size = this->paddedDim * this->cols;
		this->matrix = this->create(this->size);
		this->mask = static_cast<mask_type*>(aligned_alloc(32, this->paddedDim * sizeof(T)));

		// Offset used to accumulate instead of multiplying
		size_t dimOffset = 0;

		// registers used to peform the SIMD operations
		intrinsic_type external = Matrix2D_S<T>::intrinsic_init_register(initializer);

		// Set padding elements to 0
		for (size_t j = 0; j < this->cols; j++) {

			for (size_t i = 0; i < this->rows; i += Matrix2D_S<T>::simdStep) {
				Matrix2D_S<T>::intrinsic_store(&this->matrix[dimOffset + i], external);
			}

			for (size_t i = this->rows; i < this->paddedDim; i++) {
				this->matrix[dimOffset + i] = 0;
			}

			// Accumulate each column
			dimOffset += this->paddedDim;
		}

		// complete the mask
		for (size_t i = 0; i < this->paddedDim; i++) {
			this->mask[i] = -static_cast<mask_type>((i < this->rows));
		}

	}
}



template<typename T>
Matrix2D_S<T>::Matrix2D_S(const size_t rows, const size_t cols, const M_Orientation orientation, T(* const operation)()) {

	this->rows = rows;
	this->cols = cols;
	this->orientation = orientation;

	// The memory is allocated to the heap.
	if (orientation == M_Orientation::HORIZONTAL) {
		this->padding = Matrix2D_S<T>::simdStep - (this->cols % Matrix2D_S<T>::simdStep);
		this->padding = (this->padding == Matrix2D_S<T>::simdStep ? 0 : this->padding);
		this->paddedDim = this->cols + this->padding;
		this->size = this->paddedDim * this->rows;
		this->matrix = this->create(this->size);
		this->mask = static_cast<mask_type*>(aligned_alloc(32, this->paddedDim * sizeof(T)));

		// Offset used to accumulate instead of multiplying
		size_t dimOffset = 0;

		// Set padding elements to 0
		for (size_t i = 0; i < this->rows; i++) {

			for (size_t j = 0; j < this->cols; j++) {
				this->matrix[dimOffset + j] = operation();
			}

			for (size_t j = this->cols; j < this->paddedDim; j++) {
				this->matrix[dimOffset + j] = 0;
			}

			// Accumulate each row
			dimOffset += this->paddedDim;
		}

		// complete the mask
		for (size_t j = 0; j < this->paddedDim; j++) {
			this->mask[j] = -static_cast<mask_type>((j < this->cols));
		}

	} else if (orientation == M_Orientation::VERTICAL) {
		this->padding = Matrix2D_S<T>::simdStep - (this->rows % Matrix2D_S<T>::simdStep);
		this->padding = (this->padding == Matrix2D_S<T>::simdStep ? 0 : this->padding);
		this->paddedDim = this->rows + this->padding;
		this->size = this->paddedDim * this->cols;
		this->matrix = this->create(this->size);
		this->mask = static_cast<mask_type*>(aligned_alloc(32, this->paddedDim * sizeof(T)));

		// Offset used to accumulate instead of multiplying
		size_t dimOffset = 0;

		// Set padding elements to 0
		for (size_t j = 0; j < this->cols; j++) {

			for (size_t i = 0; i < this->rows; i++) {
				this->matrix[dimOffset + i] = operation();
			}

			for (size_t i = this->rows; i < this->paddedDim; i++) {
				this->matrix[dimOffset + i] = 0;
			}

			// Accumulate each column
			dimOffset += this->paddedDim;
		}

		// complete the mask
		for (size_t i = 0; i < this->paddedDim; i++) {
			this->mask[i] = -static_cast<mask_type>((i < this->rows));
		}

	}
}



template<typename T>
Matrix2D_S<T>::~Matrix2D_S() {
	// Delete heap allocations
	aligned_free((void*) this->matrix);
	aligned_free((void*) this->mask);
}



template<typename T>
Matrix2D_S<T> Matrix2D_S<T>::dot_H(const Matrix2D_S<T>& matrix) const {

	size_t numRows = matrix.rows;
	size_t numCols = matrix.cols;

	// Conditional throws exception for incompatible dimensions or orientations
	if (this->cols != numRows) {
		throw std::invalid_argument("The number of columns in matrix 1 does not equal the number of rows in matrix 2.");
	} else if (this->orientation != M_Orientation::HORIZONTAL || matrix.orientation != M_Orientation::VERTICAL) {
		throw std::invalid_argument("The orientation of matrix 1 must be HORIZONTAL, and the orientation of matrix 2 "
			"must be VERTICAL - for simd optimization.");
	}

	// The output matrix is created
	Matrix2D_S<T> output(this->rows, numCols, M_Orientation::HORIZONTAL);

	// registers used to peform the SIMD operations
	intrinsic_type sum;
	intrinsic_type internal;
	intrinsic_type external;

	// Offset(s) used to accumulate instead of multiplying
	size_t internalDimOffset = 0;
	size_t outputDimOffset = 0;
	size_t externalDimOffset;

	// For each row (internal)
	for (size_t i = 0; i < this->rows; i++) {

		externalDimOffset = 0;

		// Multiply and sum by all columns (external)
		for (size_t j = 0; j < numCols; j++) {

			sum = Matrix2D_S<T>::intrinsic_zeroes_register();

			// For rows (external) and columns (internal) - note strided for efficiency
			for (size_t k = 0; k < numRows; k += Matrix2D_S<T>::simdStep) {
				internal = Matrix2D_S<T>::intrinsic_load(&this->matrix[internalDimOffset + k]);
				external = Matrix2D_S<T>::intrinsic_load(&matrix.matrix[externalDimOffset + k]);
				sum = Matrix2D_S<T>::intrinsic_multiply_add(internal, external, sum);
			}

			// Store the data
			output.matrix[outputDimOffset + j] = Matrix2D_S<T>::intrinsic_aggregate_add(sum);

			// Accumulate each column (external)
			externalDimOffset += matrix.paddedDim;
		}

		// Accumulate each row internal
		internalDimOffset += this->paddedDim;
		outputDimOffset += output.paddedDim;
	}

	return output;
}



template<typename T>
Matrix2D_S<T> Matrix2D_S<T>::dot_V(const Matrix2D_S<T>& matrix) const {

	size_t numRows = matrix.rows;
	size_t numCols = matrix.cols;

	// Conditional throws exception for incompatible dimensions or orientations
	if (this->cols != numRows) {
		throw std::invalid_argument("The number of columns in matrix 1 does not equal the number of rows in matrix 2.");
	} else if (this->orientation != M_Orientation::HORIZONTAL || matrix.orientation != M_Orientation::VERTICAL) {
		throw std::invalid_argument("The orientation of matrix 1 must be HORIZONTAL, and the orientation of matrix 2 "
			"must be VERTICAL - for simd optimization.");
	}

	// The output matrix is created
	Matrix2D_S<T> output(this->rows, numCols, M_Orientation::VERTICAL);

	// registers used to peform the SIMD operations
	intrinsic_type sum;
	intrinsic_type internal;
	intrinsic_type external;

	// Offset(s) used to accumulate instead of multiplying
	size_t internalDimOffset = 0;
	size_t outputDimOffset;
	size_t externalDimOffset;

	// For each row (internal)
	for (size_t i = 0; i < this->rows; i++) {

		externalDimOffset = 0;
		outputDimOffset = 0;

		// Multiply and sum by all columns (external)
		for (size_t j = 0; j < numCols; j++) {

			sum = Matrix2D_S<T>::intrinsic_zeroes_register();

			// For rows (external) and columns (internal) - note strided for efficiency
			for (size_t k = 0; k < numRows; k += Matrix2D_S<T>::simdStep) {
				internal = Matrix2D_S<T>::intrinsic_load(&this->matrix[internalDimOffset + k]);
				external = Matrix2D_S<T>::intrinsic_load(&matrix.matrix[externalDimOffset + k]);
				sum = Matrix2D_S<T>::intrinsic_multiply_add(internal, external, sum);
			}

			// Store the data
			output.matrix[outputDimOffset + i] = Matrix2D_S<T>::intrinsic_aggregate_add(sum);

			// Accumulate each column (external)
			externalDimOffset += matrix.paddedDim;
			outputDimOffset += output.paddedDim;
		}

		// Accumulate each row internal
		internalDimOffset += this->paddedDim;
	}

	return output;
}



template<typename T>
Matrix2D_S<T> Matrix2D_S<T>::multiply(const T scalar) const {

	// The output matrix is created
	Matrix2D_S<T> output(this->rows, this->cols, this->orientation);

	// registers used to peform the SIMD operations
	intrinsic_type internal;
	intrinsic_type external = Matrix2D_S<T>::intrinsic_init_register(scalar);

	// Multiply each element by the scalar
	for (size_t i = 0; i < this->size; i += Matrix2D_S<T>::simdStep) {
		internal = Matrix2D_S<T>::intrinsic_load(&this->matrix[i]);
		Matrix2D_S<T>::intrinsic_store(&output.matrix[i], Matrix2D_S<T>::intrinsic_multiply(internal, external));
	}

	return output;
}



template<typename T>
Matrix2D_S<T> Matrix2D_S<T>::multiply(const Matrix2D_S<T>& matrix) const {

	// Validate input matrix
	if (matrix.getCols() != this->cols || matrix.getRows() != this->rows || matrix.getOrientation() != this->orientation) {
		throw std::invalid_argument("The input matrix is incompatible either in it's dimensions or orientation.");
	}

	// The output matrix is created
	Matrix2D_S<T> output(this->rows, this->cols, this->orientation);

	// registers used to peform the SIMD operations
	intrinsic_type internal;
	intrinsic_type external;

	// Multiply each element by the scalar
	for (size_t i = 0; i < this->size; i += Matrix2D_S<T>::simdStep) {
		internal = Matrix2D_S<T>::intrinsic_load(&this->matrix[i]);
		external = Matrix2D_S<T>::intrinsic_load(&matrix.matrix[i]);
		Matrix2D_S<T>::intrinsic_store(&output.matrix[i], Matrix2D_S<T>::intrinsic_multiply(internal, external));
	}

	return output;
}



template<typename T>
Matrix2D_S<T> Matrix2D_S<T>::divide(const T scalar) const {

	// The output matrix is created
	Matrix2D_S<T> output(this->rows, this->cols, this->orientation);

	// calculate the reciprocal
	T reciprocal = 1 / scalar;

	// registers used to peform the SIMD operations
	intrinsic_type internal;
	intrinsic_type external = Matrix2D_S<T>::intrinsic_init_register(reciprocal);

	// Multiply each element by the scalar
	for (size_t i = 0; i < this->size; i += Matrix2D_S<T>::simdStep) {
		internal = Matrix2D_S<T>::intrinsic_load(&this->matrix[i]);
		Matrix2D_S<T>::intrinsic_store(&output.matrix[i], Matrix2D_S<T>::intrinsic_multiply(internal, external));
	}

	return output;
}



template<typename T>
Matrix2D_S<T> Matrix2D_S<T>::divide(const Matrix2D_S<T>& matrix) const {

	// Validate input matrix
	if (matrix.getCols() != this->cols || matrix.getRows() != this->rows || matrix.getOrientation() != this->orientation) {
		throw std::invalid_argument("The input matrix is incompatible either in it's dimensions or orientation.");
	}

	// The output matrix is created
	Matrix2D_S<T> output(this->rows, this->cols, this->orientation);

	// registers used to peform the SIMD operations
	intrinsic_type internal;
	intrinsic_type external;

	// Divide each element by the corresponding element in the input
	for (size_t i = 0; i < this->size; i += Matrix2D_S<T>::simdStep) {
		internal = Matrix2D_S<T>::intrinsic_load(&this->matrix[i]);
		external = Matrix2D_S<T>::intrinsic_load(&matrix.matrix[i]);
		Matrix2D_S<T>::intrinsic_store(&output.matrix[i], Matrix2D_S<T>::intrinsic_divide(internal, external));
	}

	return output;
}



template<typename T>
Matrix2D_S<T> Matrix2D_S<T>::add(const T scalar) const {

	// The output matrix is created
	Matrix2D_S<T> output(this->rows, this->cols, this->orientation);

	// registers used to peform the SIMD operations
	intrinsic_type internal;
	intrinsic_type external = Matrix2D_S<T>::intrinsic_init_register(scalar);

	// The scalar is added to each element in the matrix
	for(size_t i = 0; i < this->size; i += Matrix2D_S<T>::simdStep) {
		internal = Matrix2D_S<T>::intrinsic_load(&this->matrix[i]);
		Matrix2D_S<T>::intrinsic_store(&output.matrix[i], Matrix2D_S<T>::intrinsic_add(internal, external));
	}

	return output;
}



template<typename T>
Matrix2D_S<T> Matrix2D_S<T>::add(const Matrix2D_S<T>& matrix) const {

	// Validate input matrix
	if (matrix.getCols() != this->cols || matrix.getRows() != this->rows || matrix.getOrientation() != this->orientation) {
		throw std::invalid_argument("The input matrix is incompatible either in it's dimensions or orientation.");
	}

	// The output matrix is created
	Matrix2D_S<T> output(this->rows, this->cols, this->orientation);

	// registers used to peform the SIMD operations
	intrinsic_type internal;
	intrinsic_type external;

	// Add each element to the corresponding element in the input matrix
	for (size_t i = 0; i < this->size; i += Matrix2D_S<T>::simdStep) {
		internal = Matrix2D_S<T>::intrinsic_load(&this->matrix[i]);
		external = Matrix2D_S<T>::intrinsic_load(&matrix.matrix[i]);
		Matrix2D_S<T>::intrinsic_store(&output.matrix[i], Matrix2D_S<T>::intrinsic_add(internal, external));
	}

	return output;
}



template<typename T>
Matrix2D_S<T> Matrix2D_S<T>::subtract(const T scalar) const {

	// The output matrix is created
	Matrix2D_S<T> output(this->rows, this->cols, this->orientation);

	// registers used to peform the SIMD operations
	intrinsic_type internal;
	intrinsic_type external = Matrix2D_S<T>::intrinsic_init_register(scalar);

	// The scalar is subtracted from every element in the matrix
	for(size_t i = 0; i < this->size; i++) {
		internal = Matrix2D_S<T>::intrinsic_load(&this->matrix[i]);
		Matrix2D_S<T>::intrinsic_store(&output.matrix[i], Matrix2D_S<T>::intrinsic_subtract(internal, external));
	}

	return output;
}



template<typename T>
Matrix2D_S<T> Matrix2D_S<T>::subtract(const Matrix2D_S<T>& matrix) const {

	// Validate input matrix
	if (matrix.getCols() != this->cols || matrix.getRows() != this->rows || matrix.getOrientation() != this->orientation) {
		throw std::invalid_argument("The input matrix is incompatible either in it's dimensions or orientation.");
	}

	// The output matrix is created
	Matrix2D_S<T> output(this->rows, this->cols, this->orientation);

	// registers used to peform the SIMD operations
	intrinsic_type internal;
	intrinsic_type external;

	// Subtract the corresponding element in the input matrix from a given element
	for (size_t i = 0; i < this->size; i += Matrix2D_S<T>::simdStep) {
		internal = Matrix2D_S<T>::intrinsic_load(&this->matrix[i]);
		external = Matrix2D_S<T>::intrinsic_load(&matrix.matrix[i]);
		Matrix2D_S<T>::intrinsic_store(&output.matrix[i], Matrix2D_S<T>::intrinsic_subtract(internal, external));
	}

	return output;
}



template<typename T>
Matrix2D_S<T> Matrix2D_S<T>::min(const T scalar) const {

	// The output matrix is created
	Matrix2D_S<T> output(this->rows, this->cols, this->orientation);

	// registers used to peform the SIMD operations
	intrinsic_type internal;
	intrinsic_type external = Matrix2D_S<T>::intrinsic_init_register(scalar);
	intrinsic_type mask;

	// Offset(s) used to accumulate instead of multiplying
	size_t dimOffset = 0;

	// The minimum is found relative to the scalar
	if (this->orientation == M_Orientation::HORIZONTAL) {

		for (size_t i = 0; i < this->rows; i++) {
			for (size_t j = 0; j < this->cols; j += Matrix2D_S<T>::simdStep) {
				mask = Matrix2D_S<T>::intrinsic_mask_load(&this->mask[j]);
				internal = Matrix2D_S<T>::intrinsic_load(&this->matrix[dimOffset + j]);
				Matrix2D_S<T>::intrinsic_store(&output.matrix[dimOffset + j], Matrix2D_S<T>::intrinsic_and(Matrix2D_S<T>::intrinsic_minimum(internal, external), mask));
			}
			dimOffset += this->paddedDim;
		}

	} else {

		for (size_t j = 0; j < this->cols; j++) {
			for (size_t i = 0; i < this->rows; i += Matrix2D_S<T>::simdStep) {
				mask = Matrix2D_S<T>::intrinsic_mask_load(&this->mask[i]);
				internal = Matrix2D_S<T>::intrinsic_load(&this->matrix[dimOffset + i]);
				Matrix2D_S<T>::intrinsic_store(&output.matrix[dimOffset + i], Matrix2D_S<T>::intrinsic_and(Matrix2D_S<T>::intrinsic_minimum(internal, external), mask));
			}
			dimOffset += this->paddedDim;
		}

	}

	return output;
}



template<typename T>
Matrix2D_S<T> Matrix2D_S<T>::min(const Matrix2D_S<T>& matrix) const {

	// Validate input matrix
	if (matrix.getCols() != this->cols || matrix.getRows() != this->rows || matrix.getOrientation() != this->orientation) {
		throw std::invalid_argument("The input matrix is incompatible either in it's dimensions or orientation.");
	}

	// The output matrix is created
	Matrix2D_S<T> output(this->rows, this->cols, this->orientation);

	// registers used to peform the SIMD operations
	intrinsic_type internal;
	intrinsic_type external;

	// Find the elementwise minimum of two matrices 
	for (size_t i = 0; i < this->size; i += Matrix2D_S<T>::simdStep) {
		internal = Matrix2D_S<T>::intrinsic_load(&this->matrix[i]);
		external = Matrix2D_S<T>::intrinsic_load(&matrix.matrix[i]);
		Matrix2D_S<T>::intrinsic_store(&output.matrix[i], Matrix2D_S<T>::intrinsic_minimum(internal, external));
	}

	return output;
}



template<typename T>
Matrix2D_S<T> Matrix2D_S<T>::max(const T scalar) const {

	// The output matrix is created
	Matrix2D_S<T> output(this->rows, this->cols, this->orientation);

	// registers used to peform the SIMD operations
	intrinsic_type internal;
	intrinsic_type external = Matrix2D_S<T>::intrinsic_init_register(scalar);
	intrinsic_type mask;

	// Offset(s) used to accumulate instead of multiplying
	size_t dimOffset = 0;

	// The maximum value is found relative to the scalar
	if (this->orientation == M_Orientation::HORIZONTAL) {

		for (size_t i = 0; i < this->rows; i++) {
			for (size_t j = 0; j < this->cols; j += Matrix2D_S<T>::simdStep) {
				mask = Matrix2D_S<T>::intrinsic_mask_load(&this->mask[j]);
				internal = Matrix2D_S<T>::intrinsic_load(&this->matrix[dimOffset + j]);
				Matrix2D_S<T>::intrinsic_store(&output.matrix[dimOffset + j], Matrix2D_S<T>::intrinsic_and(Matrix2D_S<T>::intrinsic_maximum(internal, external), mask));
			}
			dimOffset += this->paddedDim;
		}

	} else {

		for (size_t j = 0; j < this->cols; j++) {
			for (size_t i = 0; i < this->rows; i += Matrix2D_S<T>::simdStep) {
				mask = Matrix2D_S<T>::intrinsic_mask_load(&this->mask[i]);
				internal = Matrix2D_S<T>::intrinsic_load(&this->matrix[dimOffset + i]);
				Matrix2D_S<T>::intrinsic_store(&output.matrix[dimOffset + i], Matrix2D_S<T>::intrinsic_and(Matrix2D_S<T>::intrinsic_maximum(internal, external), mask));
			}
			dimOffset += this->paddedDim;
		}

	}

	return output;
}



template<typename T>
Matrix2D_S<T> Matrix2D_S<T>::max(const Matrix2D_S<T>& matrix) const {

	// Validate input matrix
	if (matrix.getCols() != this->cols || matrix.getRows() != this->rows || matrix.getOrientation() != this->orientation) {
		throw std::invalid_argument("The input matrix is incompatible either in it's dimensions or orientation.");
	}

	// The output matrix is created
	Matrix2D_S<T> output(this->rows, this->cols, this->orientation);

	// registers used to peform the SIMD operations
	intrinsic_type internal;
	intrinsic_type external;

	// Find the elementwise maximum of two matrices 
	for (size_t i = 0; i < this->size; i += Matrix2D_S<T>::simdStep) {
		internal = Matrix2D_S<T>::intrinsic_load(&this->matrix[i]);
		external = Matrix2D_S<T>::intrinsic_load(&matrix.matrix[i]);
		Matrix2D_S<T>::intrinsic_store(&output.matrix[i], Matrix2D_S<T>::intrinsic_maximum(internal, external));
	}

	return output;
}



template<typename T>
Matrix2D_S<T> Matrix2D_S<T>::apply(T(* const operation)(T)) const {

	// The output matrix is created
	Matrix2D_S<T> output(this->rows, this->cols, this->orientation);

	// The operation is performed
	for(size_t i = 0; i < this->size; i++) {
		output.matrix[i] = operation(this->matrix[i]);
	}

	return output;
}



template<typename T>
Matrix2D_S<T> Matrix2D_S<T>::transpose() const {

	// The output matrix is created
	Matrix2D_S<T> output(this->cols, this->rows, this->orientation);

	// Offset(s) used to accumulate instead of multiplying
	size_t dimOffsetInternal = 0;
	size_t dimOffsetExternal;

	// The matrix is transposed
	if(this->orientation == M_Orientation::HORIZONTAL) {
		for(size_t i = 0; i < this->rows; i++) {
			dimOffsetExternal = 0;
			for(size_t j = 0; j < this->cols; j++) {
				output.matrix[dimOffsetExternal + i] = this->matrix[dimOffsetInternal + j];
				dimOffsetExternal += output.paddedDim;
			}
			dimOffsetInternal += this->paddedDim;
		}
	} else {
		for(size_t j = 0; j < this->cols; j++) {
			dimOffsetExternal = 0;
			for(size_t i = 0; i < this->rows; i++) {
				output.matrix[dimOffsetExternal + j] = this->matrix[dimOffsetInternal + i];
				dimOffsetExternal += output.paddedDim;
			}
			dimOffsetInternal += this->paddedDim;
		}
	}

	return output;
}



template<typename T>
Matrix2D_S<T> Matrix2D_S<T>::realign() const {

	if (this->orientation == M_Orientation::HORIZONTAL) {

		Matrix2D_S<T> output(this->rows, this->cols, M_Orientation::VERTICAL);

		// Offset(s) used to accumulate instead of multiplying
		size_t dimOffsetInternal = 0;
		size_t dimOffsetExternal;

		// re-align the memory
		for (size_t i = 0; i < this->rows; i++) {
			dimOffsetExternal = 0;
			for (size_t j = 0; j < this->cols; j++) {
				output.matrix[dimOffsetExternal + i] = this->matrix[dimOffsetInternal + j];
				dimOffsetExternal += output.paddedDim;
			}
			dimOffsetInternal += this->paddedDim;
		}

		return output;

	} else {

		Matrix2D_S<T> output(this->rows, this->cols, M_Orientation::HORIZONTAL);

		// Offset(s) used to accumulate instead of multiplying
		size_t dimOffsetInternal = 0;
		size_t dimOffsetExternal;

		// re-align the memory
		for (size_t j = 0; j < this->cols; j++) {
			dimOffsetExternal = 0;
			for (size_t i = 0; i < this->rows; i++) {
				output.matrix[dimOffsetExternal + j] = this->matrix[dimOffsetInternal + i];
				dimOffsetExternal += output.paddedDim;
			}
			dimOffsetInternal += this->paddedDim;
		}

		return output;
	}
}



template<typename T>
T& Matrix2D_S<T>::operator()(const size_t row, const size_t col){

	// Validate the input row and column index
	if(row >= this->rows){
		throw std::out_of_range(std::format("Out of range row index {}.", row));
	} else if(col >= this->cols){
		throw std::out_of_range(std::format("Out of range column index {}.", col));
	}

	// Retrieve the value from the flattened array
	if(this->orientation == M_Orientation::HORIZONTAL){
		return this->matrix[row * this->paddedDim + col];
	} else{
		return this->matrix[col * this->paddedDim + row];
	}
}



template<typename T>
T& Matrix2D_S<T>::operator()(const size_t row, const size_t col) const {

	// Validate the input row and column index
	if(row >= this->rows){
		throw std::out_of_range(std::format("Out of range row index {}.", row));
	} else if(col >= this->cols){
		throw std::out_of_range(std::format("Out of range column index {}.", col));
	}

	// Retrieve the value from the flattened array
	if(this->orientation == M_Orientation::HORIZONTAL){
		return this->matrix[row * this->paddedDim + col];
	} else{
		return this->matrix[col * this->paddedDim + row];
	}
}



template<typename T>
size_t Matrix2D_S<T>::getRows() const noexcept {
	return this->rows;
}



template<typename T>
size_t Matrix2D_S<T>::getCols() const noexcept{
	return this->cols;
}



template<typename T>
size_t Matrix2D_S<T>::getPadding() const noexcept{
	return this->padding;
}



template<typename T>
M_Orientation Matrix2D_S<T>::getOrientation() const noexcept{
	return this->orientation;
}




template<typename T>
std::string Matrix2D_S<T>::toString() const noexcept{

	std::ostringstream matrix;

	if (orientation == M_Orientation::HORIZONTAL) {

		// Offset(s) used to accumulate instead of multiplying
		size_t dimOffsetInternal = 0;

		// Iterate over matrix and convert to string data.
		for (size_t i = 0; i < this->rows; i++) {
			for (size_t j = 0; j < this->cols; j++) {
				matrix << this->matrix[dimOffsetInternal + j] << ", ";
			}
			matrix << "\n";
			dimOffsetInternal += this->paddedDim;
		}

	} else {

		// Offset(s) used to accumulate instead of multiplying
		size_t dimOffsetInternal;

		// Iterate over matrix and convert to string data.
		for (size_t i = 0; i < this->rows; i++) {
			dimOffsetInternal = 0;
			for (size_t j = 0; j < this->cols; j++) {
				matrix << this->matrix[dimOffsetInternal + i] << ", ";
				dimOffsetInternal += this->paddedDim;
			}
			matrix << "\n";
		}
	}

	return matrix.str();
}