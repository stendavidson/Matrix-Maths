#include "pch.h"
#include "CppUnitTest.h"
#include "Matrix2D_S.h"

#include <vector>
#include <stdexcept>

using namespace Microsoft::VisualStudio::CppUnitTestFramework;

namespace MatrixTests
{
	TEST_CLASS(FloatMatrixTests)
	{
	public:
		
		
		TEST_METHOD(FloatMatrix_Constructor1Test_ValidVector) {

			try{
				// Set up context
				std::vector<std::vector<float>> matrixData(5, std::vector<float>(5, 0));
				Matrix2D_S<float> fiveXFive(matrixData, HORIZONTAL);

			} catch (std::invalid_argument& e) {

				std::wstring errorMessage = L"The constructor should not throw an exception for a valid 5x5 identity matrix. The following was thrown: " 
					+ std::wstring(e.what(), e.what() + std::strlen(e.what()));
				Assert::Fail(errorMessage.c_str());

			}

			// Set up context
			std::vector<std::vector<float>> matrixData(5, std::vector<float>(5, 0));
			matrixData[0][0] = 1;
			matrixData[1][1] = 1;
			matrixData[2][2] = 1;
			matrixData[3][3] = 1;
			matrixData[4][4] = 1;
			Matrix2D_S<float> fiveXFive(matrixData, HORIZONTAL);
			float** matrixPtr = fiveXFive.getMatrix();

			// Check if all values are valid
			Assert::AreEqual((long) 5, fiveXFive.getCols(), L"The number of input columns do not match");
			Assert::AreEqual((long) 5, fiveXFive.getRows(), L"The number of input columns do not match");
			Assert::AreEqual((long) 3, fiveXFive.getPadding(), L"The amount of padding does not match the correct SIMD quantity");
			Assert::AreEqual(matrixData[0][0], 1.0f, L"The matrix class has not correctly stored the input data");
			Assert::AreEqual(matrixData[1][1], 1.0f, L"The matrix class has not correctly stored the input data");
			Assert::AreEqual(matrixData[2][2], 1.0f, L"The matrix class has not correctly stored the input data");
			Assert::AreEqual(matrixData[3][3], 1.0f, L"The matrix class has not correctly stored the input data");
			Assert::AreEqual(matrixData[4][4], 1.0f, L"The matrix class has not correctly stored the input data");
		}


		TEST_METHOD(FloatMatrix_Constructor1Test_ValidVector_PaddingUnecessary) {
			
			try {
				// Set up context
				std::vector<std::vector<float>> matrixData(8, std::vector<float>(8, 0));
				Matrix2D_S<float> fiveXFive(matrixData, HORIZONTAL);

			} catch (std::invalid_argument& e) {

				std::wstring errorMessage = L"The constructor should not throw an exception for a valid 5x5 matrix. The following was thrown: " 
					+ std::wstring(e.what(), e.what() + std::strlen(e.what()));
				Assert::Fail(errorMessage.c_str());

			}

			// Set up context
			std::vector<std::vector<float>> matrixData(8, std::vector<float>(8, 0));
			Matrix2D_S<float> fiveXFive(matrixData, HORIZONTAL);

			Assert::AreEqual((long) 8, fiveXFive.getCols(), L"The number of input columns do not match");
			Assert::AreEqual((long) 8, fiveXFive.getRows(), L"The number of input columns do not match");
			Assert::AreEqual((long) 0, fiveXFive.getPadding(), L"The amount of padding does not match the correct SIMD quantity");
		}


		TEST_METHOD(FloatMatrix_Constructor1Test_InvalidVector_InconsistentRowLength) {

			try {
				// Set up context
				std::vector<std::vector<float>> matrixData = {
					std::vector<float>(8, 0),
					std::vector<float>(8, 0),
					std::vector<float>(8, 0),
					std::vector<float>(8, 0),
					std::vector<float>(8, 0),
					std::vector<float>(8, 0),
					std::vector<float>(7, 0), // Incorrect row length
					std::vector<float>(8, 0)
				};
				Matrix2D_S<float> fiveXFive(matrixData, HORIZONTAL);

			} catch (std::invalid_argument& e) {

				Assert::AreEqual(std::wstring(e.what(), e.what() + std::strlen(e.what())), std::wstring(L"Dimension inconsistencies were found in the input data."));

			}
		}


		TEST_METHOD(FloatMatrix_Constructor2Test_Valid) {

			try {

				Matrix2D_S<float> fiveXFive(5, 5, HORIZONTAL);

			} catch (std::invalid_argument& e) {

				std::wstring errorMessage = L"The constructor should not throw an exception for a valid 5x5. The following was thrown: "
					+ std::wstring(e.what(), e.what() + std::strlen(e.what()));
				Assert::Fail(errorMessage.c_str());

			}

			// Set up context
			Matrix2D_S<float> fiveXFive(5, 5, HORIZONTAL);

			Assert::AreEqual((long)5, fiveXFive.getCols(), L"The number of input columns do not match");
			Assert::AreEqual((long)5, fiveXFive.getRows(), L"The number of input columns do not match");
			Assert::AreEqual((long)3, fiveXFive.getPadding(), L"The amount of padding does not match the correct SIMD quantity");
		}

		TEST_METHOD(FloatMatrix_CopyConstructorTest_Valid) {

			try {

				Matrix2D_S<float> fiveXFive(5, 5, HORIZONTAL);
				Matrix2D_S<float> fiveXFiveCopy(fiveXFive);

			} catch (std::invalid_argument& e) {

				std::wstring errorMessage = L"The constructor should not throw an exception for a valid 5x5. The following was thrown: "
					+ std::wstring(e.what(), e.what() + std::strlen(e.what()));
				Assert::Fail(errorMessage.c_str());

			}

			// Set up context
			Matrix2D_S<float> fiveXFive(5, 5, HORIZONTAL);
			float** matrixPtr = fiveXFive.getMatrix();
			matrixPtr[0][0] = 1;
			Matrix2D_S<float> fiveXFiveCopy(fiveXFive);
			float** matrixCopyPtr = fiveXFiveCopy.getMatrix();

			Assert::AreEqual((long) 5, fiveXFive.getCols(), L"The number of input columns do not match");
			Assert::AreEqual((long) 5, fiveXFive.getRows(), L"The number of input columns do not match");
			Assert::AreEqual((long) 3, fiveXFive.getPadding(), L"The amount of padding does not match the correct SIMD quantity");
			Assert::AreEqual(matrixCopyPtr[0][0], 1.0f, L"The copy didn't correctly copy all of the data");
		}


		TEST_METHOD(FloatMatrix_Constructor3Test_Valid) {

			try {

				Matrix2D_S<float> fiveXFive(5, 5, HORIZONTAL, 4.2f);

			} catch (std::invalid_argument& e) {

				std::wstring errorMessage = L"The constructor should not throw an exception for a valid 5x5. The following was thrown: " 
					+ std::wstring(e.what(), e.what() + std::strlen(e.what()));
				Assert::Fail(errorMessage.c_str());

			}

			// Set up context
			Matrix2D_S<float> fiveXFive(5, 5, HORIZONTAL, 4.2f);
			float** matrixPtr = fiveXFive.getMatrix();

			Assert::AreEqual((long) 5, fiveXFive.getCols(), L"The number of input columns do not match");
			Assert::AreEqual((long) 5, fiveXFive.getRows(), L"The number of input columns do not match");
			Assert::AreEqual((long) 3, fiveXFive.getPadding(), L"The amount of padding does not match the correct SIMD quantity");
			Assert::AreEqual(matrixPtr[0][0], 4.2f, L"The copy didn't correctly copy all of the data");
		}


		TEST_METHOD(FloatMatrix_Constructor4Test_Valid) {

			try {

				Matrix2D_S<float> fiveXFive(5, 5, HORIZONTAL, []() { return 4.2f; });

			} catch (std::invalid_argument& e) {

				std::wstring errorMessage = L"The constructor should not throw an exception for a valid 5x5. The following was thrown: " 
					+ std::wstring(e.what(), e.what() + std::strlen(e.what()));
				Assert::Fail(errorMessage.c_str());

			}

			// Set up context
			Matrix2D_S<float> fiveXFive(5, 5, HORIZONTAL, 4.2f);
			float** matrixPtr = fiveXFive.getMatrix();

			Assert::AreEqual((long)5, fiveXFive.getCols(), L"The number of input columns do not match");
			Assert::AreEqual((long)5, fiveXFive.getRows(), L"The number of input columns do not match");
			Assert::AreEqual((long)3, fiveXFive.getPadding(), L"The amount of padding does not match the correct SIMD quantity");
			Assert::AreEqual(matrixPtr[0][0], 4.2f, L"The copy didn't correctly copy all of the data");
		}
	    

		TEST_METHOD(FloatMatrix_DotHTest_InvalidColumnsAndRows) {

			// Set up context
			Matrix2D_S<float> fiveXFive_1(5, 5, HORIZONTAL, 2.1f);
			Matrix2D_S<float> fiveXFive_2(4, 5, VERTICAL, 0.0f);
			float** matrixPtr = fiveXFive_2.getMatrix();

			// Assert No Exception
			try {

				Matrix2D_S<float> output = fiveXFive_1.dot_H(fiveXFive_2);

			} catch (std::invalid_argument& e) {

				Assert::AreEqual(std::wstring(e.what(), e.what() + std::strlen(e.what())), std::wstring(L"The number of columns in matrix 1 does not equal the"
					" number of rows in matrix 2."));

			}
		}
		

		
		TEST_METHOD(FloatMatrix_DotHTest_InvalidMemoryAlignment) {

			// Set up context
			Matrix2D_S<float> fiveXFive_1(5, 5, HORIZONTAL, 2.1f);
			Matrix2D_S<float> fiveXFive_2(5, 5, HORIZONTAL, 0.0f);
			float** matrixPtr = fiveXFive_2.getMatrix();

			// Assert No Exception
			try {

				Matrix2D_S<float> output = fiveXFive_1.dot_H(fiveXFive_2);

			} catch (std::invalid_argument& e) {

				Assert::AreEqual(std::wstring(e.what(), e.what() + std::strlen(e.what())), std::wstring(L"The orientation of matrix 1 must be horizontal, "
					L"and the orientation of matrix 2 must be vertical - for simd optimization."));

			}
		}
		
		
		TEST_METHOD(FloatMatrix_DotHTest_Valid5x5dot5x5) {

			// Set up context
			Matrix2D_S<float> fiveXFive_1(5, 5, HORIZONTAL, 2.1f);
			Matrix2D_S<float> fiveXFive_2(5, 5, VERTICAL, 0.0f);
			float** matrixPtr = fiveXFive_2.getMatrix();

			matrixPtr[0][0] = 2;
			matrixPtr[1][1] = 2;
			matrixPtr[2][2] = 2;
			matrixPtr[3][3] = 2;
			matrixPtr[4][4] = 2;

			// Assert No Exception
			try {

				Matrix2D_S<float> output = fiveXFive_1.dot_H(fiveXFive_2);
				float** outputMatrixPtr = output.getMatrix();

				Assert::AreEqual((long)5, output.getCols(), L"The number of input columns do not match");
				Assert::AreEqual((long)5, output.getRows(), L"The number of input columns do not match");
				Assert::AreEqual((long)3, output.getPadding(), L"The amount of padding does not match the correct SIMD quantity");

				for (long i = 0; i < output.getRows(); i++) {
					for (long j = 0; j < output.getCols(); j++) {
						Assert::AreEqual(4.2f, outputMatrixPtr[i][j], L"The dot product was not calculate correctly");
					}
				}

			} catch (std::invalid_argument& e) {

				std::wstring errorMessage = L"The dot_H method should not throw an exception for two correctly aligned 5x5"
					" matrices. The following was thrown: " + std::wstring(e.what(), e.what() + std::strlen(e.what()));
				Assert::Fail(errorMessage.c_str());

			}
		}

		TEST_METHOD(FloatMatrix_DotHTest_Valid5x5dot5x1) {

			// Set up context
			Matrix2D_S<float> fiveXFive_1(5, 5, HORIZONTAL, 2.1f);
			Matrix2D_S<float> fiveXFive_2(5, 1, VERTICAL, 1.0f);
			float** matrixPtr = fiveXFive_2.getMatrix();


			// Assert No Exception
			try {

				Matrix2D_S<float> output = fiveXFive_1.dot_H(fiveXFive_2);
				float** outputMatrixPtr = output.getMatrix();

				Assert::AreEqual((long) 1, output.getCols(), L"The number of input columns do not match");
				Assert::AreEqual((long) 5, output.getRows(), L"The number of input columns do not match");
				Assert::AreEqual((long) 7, output.getPadding(), L"The amount of padding does not match the correct SIMD quantity");

				for (long i = 0; i < output.getRows(); i++) {
					for (long j = 0; j < output.getCols(); j++) {
						Assert::AreEqual(10.5f, outputMatrixPtr[i][j], L"The dot product was not calculate correctly");
					}
				}

			} catch (std::invalid_argument& e) {

				std::wstring errorMessage = L"The dot_H method should not throw an exception for two correctly aligned 5x5"
					" matrices. The following was thrown: " + std::wstring(e.what(), e.what() + std::strlen(e.what()));
				Assert::Fail(errorMessage.c_str());

			}
		}


		TEST_METHOD(FloatMatrix_DotVTest_Valid5x5dot5x5) {

			// Set up context
			Matrix2D_S<float> fiveXFive_1(5, 5, HORIZONTAL, 2.1f);
			Matrix2D_S<float> fiveXFive_2(5, 5, VERTICAL, 0.0f);
			float** matrixPtr = fiveXFive_2.getMatrix();

			matrixPtr[0][0] = 2;
			matrixPtr[1][1] = 2;
			matrixPtr[2][2] = 2;
			matrixPtr[3][3] = 2;
			matrixPtr[4][4] = 2;

			// Assert No Exception
			try {

				Matrix2D_S<float> output = fiveXFive_1.dot_V(fiveXFive_2);
				float** outputMatrixPtr = output.getMatrix();

				Assert::AreEqual((long)5, output.getCols(), L"The number of input columns do not match");
				Assert::AreEqual((long)5, output.getRows(), L"The number of input columns do not match");
				Assert::AreEqual((long)3, output.getPadding(), L"The amount of padding does not match the correct SIMD quantity");

				for (long i = 0; i < output.getRows(); i++) {
					for (long j = 0; j < output.getCols(); j++) {
						Assert::AreEqual(4.2f, outputMatrixPtr[j][i], L"The dot product was not calculate correctly");
					}
				}

			} catch (std::invalid_argument& e) {

				std::wstring errorMessage = L"The dot_V method should not throw an exception for two correctly aligned 5x5"
					" matrices. The following was thrown: " + std::wstring(e.what(), e.what() + std::strlen(e.what()));
				Assert::Fail(errorMessage.c_str());

			}
		}



		TEST_METHOD(FloatMatrix_DotVTest_Valid5x5dot5x1) {

			// Set up context
			Matrix2D_S<float> fiveXFive_1(5, 5, HORIZONTAL, 2.1f);
			Matrix2D_S<float> fiveXFive_2(5, 1, VERTICAL, 1.0f);
			float** matrixPtr = fiveXFive_2.getMatrix();


			// Assert No Exception
			try {

				Matrix2D_S<float> output = fiveXFive_1.dot_V(fiveXFive_2);
				float** outputMatrixPtr = output.getMatrix();

				Assert::AreEqual((long) 1, output.getCols(), L"The number of input columns do not match");
				Assert::AreEqual((long) 5, output.getRows(), L"The number of input columns do not match");
				Assert::AreEqual((long) 3, output.getPadding(), L"The amount of padding does not match the correct SIMD quantity");

				for (long i = 0; i < output.getRows(); i++) {
					for (long j = 0; j < output.getCols(); j++) {
						Assert::AreEqual(10.5f, outputMatrixPtr[j][i], L"The dot product was not calculate correctly");
					}
				}

			} catch (std::invalid_argument& e) {

				std::wstring errorMessage = L"The dot_V method should not throw an exception for two correctly aligned 5x5"
					" matrices. The following was thrown: " + std::wstring(e.what(), e.what() + std::strlen(e.what()));
				Assert::Fail(errorMessage.c_str());

			}
		}


		TEST_METHOD(FloatMatrix_MultiplyTest_Valid5x5xScalar) {

			// Set up context
			Matrix2D_S<float> fiveXFive_1(5, 5, HORIZONTAL, 2.0f);

			// Assert No Exception
			try {

				Matrix2D_S<float> output = fiveXFive_1.multiply(3.0f);
				float** outputMatrixPtr = output.getMatrix();

				Assert::AreEqual((long) 5, output.getCols(), L"The number of input columns do not match");
				Assert::AreEqual((long) 5, output.getRows(), L"The number of input columns do not match");
				Assert::AreEqual((long) 3, output.getPadding(), L"The amount of padding does not match the correct SIMD quantity");

				for (long i = 0; i < output.getRows(); i++) {
					for (long j = 0; j < output.getCols(); j++) {
						Assert::AreEqual(6, (int) outputMatrixPtr[i][j], L"The multiplication was not calculate correctly");
					}
				}

			} catch (std::invalid_argument& e) {

				std::wstring errorMessage = L"The multiply method should not throw an exception for a 5x5 matrices multiplied"
					" by a scalar. The following was thrown: " + std::wstring(e.what(), e.what() + std::strlen(e.what()));
				Assert::Fail(errorMessage.c_str());

			}
		}


		TEST_METHOD(FloatMatrix_MultiplyTest_Valid5x5x5x5) {

			// Set up context
			Matrix2D_S<float> fiveXFive_1(5, 5, HORIZONTAL, 2.0f);
			Matrix2D_S<float> fiveXFive_2(5, 5, HORIZONTAL, 2.0f);

			// Assert No Exception
			try {

				Matrix2D_S<float> output = fiveXFive_1.multiply(fiveXFive_2);
				float** outputMatrixPtr = output.getMatrix();

				Assert::AreEqual((long) 5, output.getCols(), L"The number of input columns do not match");
				Assert::AreEqual((long) 5, output.getRows(), L"The number of input columns do not match");
				Assert::AreEqual((long) 3, output.getPadding(), L"The amount of padding does not match the correct SIMD quantity");

				for (long i = 0; i < output.getRows(); i++) {
					for (long j = 0; j < output.getCols(); j++) {
						Assert::AreEqual(4, (int) outputMatrixPtr[i][j], L"The multiplication was not calculate correctly");
					}
				}

			} catch (std::invalid_argument& e) {

				std::wstring errorMessage = L"The multiply method should not throw an exception for a 5x5 matrices multiplied"
					" by a matrix. The following was thrown: " + std::wstring(e.what(), e.what() + std::strlen(e.what()));
				Assert::Fail(errorMessage.c_str());

			}
		}


		TEST_METHOD(FloatMatrix_DivideTest_Valid5x5xScalar) {

			// Set up context
			Matrix2D_S<float> fiveXFive_1(5, 5, HORIZONTAL, 2.0f);

			// Assert No Exception
			try {

				Matrix2D_S<float> output = fiveXFive_1.divide(2.0f);
				float** outputMatrixPtr = output.getMatrix();

				Assert::AreEqual((long) 5, output.getCols(), L"The number of input columns do not match");
				Assert::AreEqual((long) 5, output.getRows(), L"The number of input columns do not match");
				Assert::AreEqual((long) 3, output.getPadding(), L"The amount of padding does not match the correct SIMD quantity");

				for (long i = 0; i < output.getRows(); i++) {
					for (long j = 0; j < output.getCols(); j++) {
						Assert::AreEqual(1, (int)outputMatrixPtr[i][j], L"The division was not calculate correctly");
					}
				}

			} catch (std::invalid_argument& e) {

				std::wstring errorMessage = L"The division method should not throw an exception for a 5x5 matrices divided"
					" by a scalar. The following was thrown: " + std::wstring(e.what(), e.what() + std::strlen(e.what()));
				Assert::Fail(errorMessage.c_str());

			}
		}


		TEST_METHOD(FloatMatrix_DivideTest_Valid5x5x5x5) {

			// Set up context
			Matrix2D_S<float> fiveXFive_1(5, 5, HORIZONTAL, 2.0f);
			Matrix2D_S<float> fiveXFive_2(5, 5, HORIZONTAL, 2.0f);

			// Assert No Exception
			try {

				Matrix2D_S<float> output = fiveXFive_1.divide(fiveXFive_2);
				float** outputMatrixPtr = output.getMatrix();

				Assert::AreEqual((long) 5, output.getCols(), L"The number of input columns do not match");
				Assert::AreEqual((long) 5, output.getRows(), L"The number of input columns do not match");
				Assert::AreEqual((long) 3, output.getPadding(), L"The amount of padding does not match the correct SIMD quantity");

				for (long i = 0; i < output.getRows(); i++) {
					for (long j = 0; j < output.getCols(); j++) {
						Assert::AreEqual(1, (int) outputMatrixPtr[i][j], L"The division was not calculate correctly");
					}
				}

			} catch (std::invalid_argument& e) {

				std::wstring errorMessage = L"The divide method should not throw an exception for a 5x5 matrices divided"
					" by a matrix. The following was thrown: " + std::wstring(e.what(), e.what() + std::strlen(e.what()));
				Assert::Fail(errorMessage.c_str());

			}
		}


		TEST_METHOD(FloatMatrix_AddTest_Valid5x5xScalar) {

			// Set up context
			Matrix2D_S<float> fiveXFive_1(5, 5, HORIZONTAL, 2.0f);

			// Assert No Exception
			try {

				Matrix2D_S<float> output = fiveXFive_1.add(2.0f);
				float** outputMatrixPtr = output.getMatrix();

				Assert::AreEqual((long) 5, output.getCols(), L"The number of input columns do not match");
				Assert::AreEqual((long) 5, output.getRows(), L"The number of input columns do not match");
				Assert::AreEqual((long) 3, output.getPadding(), L"The amount of padding does not match the correct SIMD quantity");

				for (long i = 0; i < output.getRows(); i++) {
					for (long j = 0; j < output.getCols(); j++) {
						Assert::AreEqual(4, (int) outputMatrixPtr[i][j], L"The addition was not calculate correctly");
					}
				}

			} catch (std::invalid_argument& e) {

				std::wstring errorMessage = L"The add method should not throw an exception for a 5x5 matrices divided"
					" by a scalar. The following was thrown: " + std::wstring(e.what(), e.what() + std::strlen(e.what()));
				Assert::Fail(errorMessage.c_str());

			}
		}


		TEST_METHOD(FloatMatrix_AddTest_Valid5x5x5x5) {

			// Set up context
			Matrix2D_S<float> fiveXFive_1(5, 5, HORIZONTAL, 2.0f);
			Matrix2D_S<float> fiveXFive_2(5, 5, HORIZONTAL, 2.0f);

			// Assert No Exception
			try {

				Matrix2D_S<float> output = fiveXFive_1.add(fiveXFive_2);
				float** outputMatrixPtr = output.getMatrix();

				Assert::AreEqual((long) 5, output.getCols(), L"The number of input columns do not match");
				Assert::AreEqual((long) 5, output.getRows(), L"The number of input columns do not match");
				Assert::AreEqual((long) 3, output.getPadding(), L"The amount of padding does not match the correct SIMD quantity");

				for (long i = 0; i < output.getRows(); i++) {
					for (long j = 0; j < output.getCols(); j++) {
						Assert::AreEqual(4, (int) outputMatrixPtr[i][j], L"The addition was not calculate correctly");
					}
				}

			} catch (std::invalid_argument& e) {

				std::wstring errorMessage = L"The add method should not throw an exception for a 5x5 matrices divided"
					" by a matrix. The following was thrown: " + std::wstring(e.what(), e.what() + std::strlen(e.what()));
				Assert::Fail(errorMessage.c_str());

			}
		}


		TEST_METHOD(FloatMatrix_SubtractTest_Valid5x5xScalar) {

			// Set up context
			Matrix2D_S<float> fiveXFive_1(5, 5, HORIZONTAL, 4.0f);

			// Assert No Exception
			try {

				Matrix2D_S<float> output = fiveXFive_1.subtract(2.0f);
				float** outputMatrixPtr = output.getMatrix();

				Assert::AreEqual((long) 5, output.getCols(), L"The number of input columns do not match");
				Assert::AreEqual((long) 5, output.getRows(), L"The number of input columns do not match");
				Assert::AreEqual((long) 3, output.getPadding(), L"The amount of padding does not match the correct SIMD quantity");

				for (long i = 0; i < output.getRows(); i++) {
					for (long j = 0; j < output.getCols(); j++) {
						Assert::AreEqual(2, (int) outputMatrixPtr[i][j], L"The subtraction was not calculate correctly");
					}
				}

			} catch (std::invalid_argument& e) {

				std::wstring errorMessage = L"The subtract method should not throw an exception for a 5x5 matrices subtracted"
					" by a scalar. The following was thrown: " + std::wstring(e.what(), e.what() + std::strlen(e.what()));
				Assert::Fail(errorMessage.c_str());

			}
		}


		TEST_METHOD(FloatMatrix_SubtractTest_Valid5x5x5x5) {

			// Set up context
			Matrix2D_S<float> fiveXFive_1(5, 5, HORIZONTAL, 4.0f);
			Matrix2D_S<float> fiveXFive_2(5, 5, HORIZONTAL, 2.0f);

			// Assert No Exception
			try {

				Matrix2D_S<float> output = fiveXFive_1.subtract(fiveXFive_2);
				float** outputMatrixPtr = output.getMatrix();

				Assert::AreEqual((long) 5, output.getCols(), L"The number of input columns do not match");
				Assert::AreEqual((long) 5, output.getRows(), L"The number of input columns do not match");
				Assert::AreEqual((long) 3, output.getPadding(), L"The amount of padding does not match the correct SIMD quantity");

				for (long i = 0; i < output.getRows(); i++) {
					for (long j = 0; j < output.getCols(); j++) {
						Assert::AreEqual(2, (int) outputMatrixPtr[i][j], L"The subtraction was not calculate correctly");
					}
				}

			} catch (std::invalid_argument& e) {

				std::wstring errorMessage = L"The subtract method should not throw an exception for a 5x5 matrices subtracted"
					" by a matrix. The following was thrown: " + std::wstring(e.what(), e.what() + std::strlen(e.what()));
				Assert::Fail(errorMessage.c_str());

			}
		}


		TEST_METHOD(FloatMatrix_MinTest_Valid5x5xScalar) {

			// Set up context
			Matrix2D_S<float> fiveXFive_1(5, 5, HORIZONTAL, 4.0f);

			// Assert No Exception
			try {

				Matrix2D_S<float> output = fiveXFive_1.min(2.0f);
				float** outputMatrixPtr = output.getMatrix();

				Assert::AreEqual((long) 5, output.getCols(), L"The number of input columns do not match");
				Assert::AreEqual((long) 5, output.getRows(), L"The number of input columns do not match");
				Assert::AreEqual((long) 3, output.getPadding(), L"The amount of padding does not match the correct SIMD quantity");

				// Validate that the minimum value was found correctly
				for (long i = 0; i < output.getRows(); i++) {
					for (long j = 0; j < output.getCols(); j++) {
						Assert::AreEqual(2, (int) outputMatrixPtr[i][j], L"The minimum was not found correctly");
					}
				}

				for (long i = 0; i < output.getRows(); i++) {
					for (long j = output.getCols(); j < output.getCols() + output.getPadding(); j++) {
						Assert::AreEqual(0, (int) outputMatrixPtr[i][j], L"The minimum function was applied incorrectly to SIMD padding.");
					}
				}

			} catch (std::invalid_argument& e) {

				std::wstring errorMessage = L"The min method should not throw an exception for a 5x5 matrice and"
					" a scalar. The following was thrown: " + std::wstring(e.what(), e.what() + std::strlen(e.what()));
				Assert::Fail(errorMessage.c_str());

			}
		}


		TEST_METHOD(FloatMatrix_MinTest_Valid5x5x5x5) {

			// Set up context
			Matrix2D_S<float> fiveXFive_1(5, 5, HORIZONTAL, 4.0f);
			Matrix2D_S<float> fiveXFive_2(5, 5, HORIZONTAL, 2.0f);

			// Assert No Exception
			try {

				Matrix2D_S<float> output = fiveXFive_1.min(fiveXFive_2);
				float** outputMatrixPtr = output.getMatrix();

				Assert::AreEqual((long)5, output.getCols(), L"The number of input columns do not match");
				Assert::AreEqual((long)5, output.getRows(), L"The number of input columns do not match");
				Assert::AreEqual((long)3, output.getPadding(), L"The amount of padding does not match the correct SIMD quantity");

				for (long i = 0; i < output.getRows(); i++) {
					for (long j = 0; j < output.getCols(); j++) {
						Assert::AreEqual(2, (int) outputMatrixPtr[i][j], L"The minimum was not found correctly");
					}
				}

			} catch (std::invalid_argument& e) {

				std::wstring errorMessage = L"The min method should not throw an exception for two 5x5 matrices"
					". The following was thrown: " + std::wstring(e.what(), e.what() + std::strlen(e.what()));
				Assert::Fail(errorMessage.c_str());

			}
		}


		TEST_METHOD(FloatMatrix_MaxTest_Valid5x5xScalar) {

			// Set up context
			Matrix2D_S<float> fiveXFive_1(5, 5, HORIZONTAL, 0.0f);

			// Assert No Exception
			try {

				Matrix2D_S<float> output = fiveXFive_1.max(2.0f);
				float** outputMatrixPtr = output.getMatrix();

				Assert::AreEqual((long) 5, output.getCols(), L"The number of input columns do not match");
				Assert::AreEqual((long) 5, output.getRows(), L"The number of input columns do not match");
				Assert::AreEqual((long) 3, output.getPadding(), L"The amount of padding does not match the correct SIMD quantity");

				// Validate that the minimum value was found correctly
				for (long i = 0; i < output.getRows(); i++) {
					for (long j = 0; j < output.getCols(); j++) {
						Assert::AreEqual(2, (int) outputMatrixPtr[i][j], L"The maximum was not found correctly");
					}
				}

				for (long i = 0; i < output.getRows(); i++) {
					for (long j = output.getCols(); j < output.getCols() + output.getPadding(); j++) {
						Assert::AreEqual(0, (int) outputMatrixPtr[i][j], L"The maximum function was applied incorrectly to SIMD padding.");
					}
				}

			} catch (std::invalid_argument& e) {

				std::wstring errorMessage = L"The max method should not throw an exception for a 5x5 matrice and"
					" a scalar. The following was thrown: " + std::wstring(e.what(), e.what() + std::strlen(e.what()));
				Assert::Fail(errorMessage.c_str());

			}
		}


		TEST_METHOD(FloatMatrix_MaxTest_Valid5x5x5x5) {

			// Set up context
			Matrix2D_S<float> fiveXFive_1(5, 5, HORIZONTAL, 4.0f);
			Matrix2D_S<float> fiveXFive_2(5, 5, HORIZONTAL, 2.0f);

			// Assert No Exception
			try {

				Matrix2D_S<float> output = fiveXFive_1.max(fiveXFive_2);
				float** outputMatrixPtr = output.getMatrix();

				Assert::AreEqual((long) 5, output.getCols(), L"The number of input columns do not match");
				Assert::AreEqual((long) 5, output.getRows(), L"The number of input columns do not match");
				Assert::AreEqual((long) 3, output.getPadding(), L"The amount of padding does not match the correct SIMD quantity");

				for (long i = 0; i < output.getRows(); i++) {
					for (long j = 0; j < output.getCols(); j++) {
						Assert::AreEqual(4, (int)outputMatrixPtr[i][j], L"The maximum was not found correctly");
					}
				}

			} catch (std::invalid_argument& e) {

				std::wstring errorMessage = L"The max method should not throw an exception for two 5x5 matrices"
					". The following was thrown: " + std::wstring(e.what(), e.what() + std::strlen(e.what()));
				Assert::Fail(errorMessage.c_str());

			}
		}

		TEST_METHOD(FloatMatrixV_MultiplyTest_Valid5x5xScalar) {

			// Set up context
			Matrix2D_S<float> fiveXFive_1(5, 5, VERTICAL, 2.0f);

			// Assert No Exception
			try {

				Matrix2D_S<float> output = fiveXFive_1.multiply(3.0f);
				float** outputMatrixPtr = output.getMatrix();

				Assert::AreEqual((long)5, output.getCols(), L"The number of input columns do not match");
				Assert::AreEqual((long)5, output.getRows(), L"The number of input columns do not match");
				Assert::AreEqual((long)3, output.getPadding(), L"The amount of padding does not match the correct SIMD quantity");

				for (long i = 0; i < output.getRows(); i++) {
					for (long j = 0; j < output.getCols(); j++) {
						Assert::AreEqual(6, (int)outputMatrixPtr[i][j], L"The multiplication was not calculate correctly");
					}
				}

			} catch (std::invalid_argument& e) {

				std::wstring errorMessage = L"The multiply method should not throw an exception for a 5x5 matrices multiplied"
					" by a scalar. The following was thrown: " + std::wstring(e.what(), e.what() + std::strlen(e.what()));
				Assert::Fail(errorMessage.c_str());

			}
		}


		TEST_METHOD(FloatMatrixV_MultiplyTest_Valid5x5x5x5) {

			// Set up context
			Matrix2D_S<float> fiveXFive_1(5, 5, VERTICAL, 2.0f);
			Matrix2D_S<float> fiveXFive_2(5, 5, VERTICAL, 2.0f);

			// Assert No Exception
			try {

				Matrix2D_S<float> output = fiveXFive_1.multiply(fiveXFive_2);
				float** outputMatrixPtr = output.getMatrix();

				Assert::AreEqual((long)5, output.getCols(), L"The number of input columns do not match");
				Assert::AreEqual((long)5, output.getRows(), L"The number of input columns do not match");
				Assert::AreEqual((long)3, output.getPadding(), L"The amount of padding does not match the correct SIMD quantity");

				for (long i = 0; i < output.getRows(); i++) {
					for (long j = 0; j < output.getCols(); j++) {
						Assert::AreEqual(4, (int)outputMatrixPtr[i][j], L"The multiplication was not calculate correctly");
					}
				}

			} catch (std::invalid_argument& e) {

				std::wstring errorMessage = L"The multiply method should not throw an exception for a 5x5 matrices multiplied"
					" by a matrix. The following was thrown: " + std::wstring(e.what(), e.what() + std::strlen(e.what()));
				Assert::Fail(errorMessage.c_str());

			}
		}


		TEST_METHOD(FloatMatrixV_DivideTest_Valid5x5xScalar) {

			// Set up context
			Matrix2D_S<float> fiveXFive_1(5, 5, VERTICAL, 2.0f);

			// Assert No Exception
			try {

				Matrix2D_S<float> output = fiveXFive_1.divide(2.0f);
				float** outputMatrixPtr = output.getMatrix();

				Assert::AreEqual((long)5, output.getCols(), L"The number of input columns do not match");
				Assert::AreEqual((long)5, output.getRows(), L"The number of input columns do not match");
				Assert::AreEqual((long)3, output.getPadding(), L"The amount of padding does not match the correct SIMD quantity");

				for (long i = 0; i < output.getRows(); i++) {
					for (long j = 0; j < output.getCols(); j++) {
						Assert::AreEqual(1, (int)outputMatrixPtr[i][j], L"The division was not calculate correctly");
					}
				}

			} catch (std::invalid_argument& e) {

				std::wstring errorMessage = L"The divide method should not throw an exception for a 5x5 matrices divided"
					" by a scalar. The following was thrown: " + std::wstring(e.what(), e.what() + std::strlen(e.what()));
				Assert::Fail(errorMessage.c_str());

			}
		}


		TEST_METHOD(FloatMatrixV_DivideTest_Valid5x5x5x5) {

			// Set up context
			Matrix2D_S<float> fiveXFive_1(5, 5, VERTICAL, 2.0f);
			Matrix2D_S<float> fiveXFive_2(5, 5, VERTICAL, 2.0f);

			// Assert No Exception
			try {

				Matrix2D_S<float> output = fiveXFive_1.divide(fiveXFive_2);
				float** outputMatrixPtr = output.getMatrix();

				Assert::AreEqual((long)5, output.getCols(), L"The number of input columns do not match");
				Assert::AreEqual((long)5, output.getRows(), L"The number of input columns do not match");
				Assert::AreEqual((long)3, output.getPadding(), L"The amount of padding does not match the correct SIMD quantity");

				for (long i = 0; i < output.getRows(); i++) {
					for (long j = 0; j < output.getCols(); j++) {
						Assert::AreEqual(1, (int)outputMatrixPtr[i][j], L"The division was not calculate correctly");
					}
				}

			} catch (std::invalid_argument& e) {

				std::wstring errorMessage = L"The divide method should not throw an exception for a 5x5 matrices divided"
					" by a matrix. The following was thrown: " + std::wstring(e.what(), e.what() + std::strlen(e.what()));
				Assert::Fail(errorMessage.c_str());

			}
		}


		TEST_METHOD(FloatMatrixV_AddTest_Valid5x5xScalar) {

			// Set up context
			Matrix2D_S<float> fiveXFive_1(5, 5, VERTICAL, 2.0f);

			// Assert No Exception
			try {

				Matrix2D_S<float> output = fiveXFive_1.add(2.0f);
				float** outputMatrixPtr = output.getMatrix();

				Assert::AreEqual((long)5, output.getCols(), L"The number of input columns do not match");
				Assert::AreEqual((long)5, output.getRows(), L"The number of input columns do not match");
				Assert::AreEqual((long)3, output.getPadding(), L"The amount of padding does not match the correct SIMD quantity");

				for (long i = 0; i < output.getRows(); i++) {
					for (long j = 0; j < output.getCols(); j++) {
						Assert::AreEqual(4, (int)outputMatrixPtr[i][j], L"The addition was not calculate correctly");
					}
				}

			} catch (std::invalid_argument& e) {

				std::wstring errorMessage = L"The add method should not throw an exception for a 5x5 matrices divided"
					" by a scalar. The following was thrown: " + std::wstring(e.what(), e.what() + std::strlen(e.what()));
				Assert::Fail(errorMessage.c_str());

			}
		}


		TEST_METHOD(FloatMatrixV_AddTest_Valid5x5x5x5) {

			// Set up context
			Matrix2D_S<float> fiveXFive_1(5, 5, VERTICAL, 2.0f);
			Matrix2D_S<float> fiveXFive_2(5, 5, VERTICAL, 2.0f);

			// Assert No Exception
			try {

				Matrix2D_S<float> output = fiveXFive_1.add(fiveXFive_2);
				float** outputMatrixPtr = output.getMatrix();

				Assert::AreEqual((long)5, output.getCols(), L"The number of input columns do not match");
				Assert::AreEqual((long)5, output.getRows(), L"The number of input columns do not match");
				Assert::AreEqual((long)3, output.getPadding(), L"The amount of padding does not match the correct SIMD quantity");

				for (long i = 0; i < output.getRows(); i++) {
					for (long j = 0; j < output.getCols(); j++) {
						Assert::AreEqual(4, (int)outputMatrixPtr[i][j], L"The addition was not calculate correctly");
					}
				}

			} catch (std::invalid_argument& e) {

				std::wstring errorMessage = L"The add method should not throw an exception for a 5x5 matrices divided"
					" by a matrix. The following was thrown: " + std::wstring(e.what(), e.what() + std::strlen(e.what()));
				Assert::Fail(errorMessage.c_str());

			}
		}


		TEST_METHOD(FloatMatrixV_SubtractTest_Valid5x5xScalar) {

			// Set up context
			Matrix2D_S<float> fiveXFive_1(5, 5, VERTICAL, 4.0f);

			// Assert No Exception
			try {

				Matrix2D_S<float> output = fiveXFive_1.subtract(2.0f);
				float** outputMatrixPtr = output.getMatrix();

				Assert::AreEqual((long)5, output.getCols(), L"The number of input columns do not match");
				Assert::AreEqual((long)5, output.getRows(), L"The number of input columns do not match");
				Assert::AreEqual((long)3, output.getPadding(), L"The amount of padding does not match the correct SIMD quantity");

				for (long i = 0; i < output.getRows(); i++) {
					for (long j = 0; j < output.getCols(); j++) {
						Assert::AreEqual(2, (int)outputMatrixPtr[i][j], L"The subtraction was not calculate correctly");
					}
				}

			} catch (std::invalid_argument& e) {

				std::wstring errorMessage = L"The subtract method should not throw an exception for a 5x5 matrices subtracted"
					" by a scalar. The following was thrown: " + std::wstring(e.what(), e.what() + std::strlen(e.what()));
				Assert::Fail(errorMessage.c_str());

			}
		}


		TEST_METHOD(FloatMatrixV_SubtractTest_Valid5x5x5x5) {

			// Set up context
			Matrix2D_S<float> fiveXFive_1(5, 5, VERTICAL, 4.0f);
			Matrix2D_S<float> fiveXFive_2(5, 5, VERTICAL, 2.0f);

			// Assert No Exception
			try {

				Matrix2D_S<float> output = fiveXFive_1.subtract(fiveXFive_2);
				float** outputMatrixPtr = output.getMatrix();

				Assert::AreEqual((long)5, output.getCols(), L"The number of input columns do not match");
				Assert::AreEqual((long)5, output.getRows(), L"The number of input columns do not match");
				Assert::AreEqual((long)3, output.getPadding(), L"The amount of padding does not match the correct SIMD quantity");

				for (long i = 0; i < output.getRows(); i++) {
					for (long j = 0; j < output.getCols(); j++) {
						Assert::AreEqual(2, (int)outputMatrixPtr[i][j], L"The subtraction was not calculate correctly");
					}
				}

			} catch (std::invalid_argument& e) {

				std::wstring errorMessage = L"The subtract method should not throw an exception for a 5x5 matrices subtracted"
					" by a matrix. The following was thrown: " + std::wstring(e.what(), e.what() + std::strlen(e.what()));
				Assert::Fail(errorMessage.c_str());

			}
		}


		TEST_METHOD(FloatMatrixV_MinTest_Valid5x5xScalar) {

			// Set up context
			Matrix2D_S<float> fiveXFive_1(5, 5, VERTICAL, 4.0f);

			// Assert No Exception
			try {

				Matrix2D_S<float> output = fiveXFive_1.min(2.0f);
				float** outputMatrixPtr = output.getMatrix();

				Assert::AreEqual((long)5, output.getCols(), L"The number of input columns do not match");
				Assert::AreEqual((long)5, output.getRows(), L"The number of input columns do not match");
				Assert::AreEqual((long)3, output.getPadding(), L"The amount of padding does not match the correct SIMD quantity");

				for (long i = 0; i < output.getRows(); i++) {
					for (long j = 0; j < output.getCols(); j++) {
						Assert::AreEqual(2, (int)outputMatrixPtr[i][j], L"The minimum was not found correctly");
					}
				}

				for (long i = 0; i < output.getRows(); i++) {
					for (long j = output.getCols(); j < output.getCols() + output.getPadding(); j++) {
						Assert::AreEqual(0, (int)outputMatrixPtr[i][j], L"The minimum function was applied incorrectly to SIMD padding.");
					}
				}

			} catch (std::invalid_argument& e) {

				std::wstring errorMessage = L"The min method should not throw an exception for a 5x5 matrice and"
					" a scalar. The following was thrown: " + std::wstring(e.what(), e.what() + std::strlen(e.what()));
				Assert::Fail(errorMessage.c_str());

			}
		}


		TEST_METHOD(FloatMatrixV_MinTest_Valid5x5x5x5) {

			// Set up context
			Matrix2D_S<float> fiveXFive_1(5, 5, VERTICAL, 4.0f);
			Matrix2D_S<float> fiveXFive_2(5, 5, VERTICAL, 2.0f);

			// Assert No Exception
			try {

				Matrix2D_S<float> output = fiveXFive_1.min(fiveXFive_2);
				float** outputMatrixPtr = output.getMatrix();

				Assert::AreEqual((long)5, output.getCols(), L"The number of input columns do not match");
				Assert::AreEqual((long)5, output.getRows(), L"The number of input columns do not match");
				Assert::AreEqual((long)3, output.getPadding(), L"The amount of padding does not match the correct SIMD quantity");

				for (long i = 0; i < output.getRows(); i++) {
					for (long j = 0; j < output.getCols(); j++) {
						Assert::AreEqual(2, (int)outputMatrixPtr[i][j], L"The minimum was not found correctly");
					}
				}

			} catch (std::invalid_argument& e) {

				std::wstring errorMessage = L"The min method should not throw an exception for two 5x5 matrices"
					". The following was thrown: " + std::wstring(e.what(), e.what() + std::strlen(e.what()));
				Assert::Fail(errorMessage.c_str());

			}
		}


		TEST_METHOD(FloatMatrixV_MaxTest_Valid5x5xScalar) {

			// Set up context
			Matrix2D_S<float> fiveXFive_1(5, 5, VERTICAL, 0.0f);

			// Assert No Exception
			try {

				Matrix2D_S<float> output = fiveXFive_1.max(2.0f);
				float** outputMatrixPtr = output.getMatrix();

				Assert::AreEqual((long)5, output.getCols(), L"The number of input columns do not match");
				Assert::AreEqual((long)5, output.getRows(), L"The number of input columns do not match");
				Assert::AreEqual((long)3, output.getPadding(), L"The amount of padding does not match the correct SIMD quantity");

				for (long i = 0; i < output.getRows(); i++) {
					for (long j = 0; j < output.getCols(); j++) {
						Assert::AreEqual(2, (int)outputMatrixPtr[i][j], L"The maximum was not found correctly");
					}
				}

				for (long i = 0; i < output.getRows(); i++) {
					for (long j = output.getCols(); j < output.getCols() + output.getPadding(); j++) {
						Assert::AreEqual(0, (int)outputMatrixPtr[i][j], L"The maximum function was applied incorrectly to SIMD padding.");
					}
				}

			} catch (std::invalid_argument& e) {

				std::wstring errorMessage = L"The max method should not throw an exception for a 5x5 matrice and"
					" a scalar. The following was thrown: " + std::wstring(e.what(), e.what() + std::strlen(e.what()));
				Assert::Fail(errorMessage.c_str());

			}
		}


		TEST_METHOD(FloatMatrixV_MaxTest_Valid5x5x5x5) {

			// Set up context
			Matrix2D_S<float> fiveXFive_1(5, 5, VERTICAL, 4.0f);
			Matrix2D_S<float> fiveXFive_2(5, 5, VERTICAL, 2.0f);

			// Assert No Exception
			try {

				Matrix2D_S<float> output = fiveXFive_1.max(fiveXFive_2);
				float** outputMatrixPtr = output.getMatrix();

				Assert::AreEqual((long)5, output.getCols(), L"The number of input columns do not match");
				Assert::AreEqual((long)5, output.getRows(), L"The number of input columns do not match");
				Assert::AreEqual((long)3, output.getPadding(), L"The amount of padding does not match the correct SIMD quantity");

				for (long i = 0; i < output.getRows(); i++) {
					for (long j = 0; j < output.getCols(); j++) {
						Assert::AreEqual(4, (int)outputMatrixPtr[i][j], L"The maximum was not found correctly");
					}
				}

			} catch (std::invalid_argument& e) {

				std::wstring errorMessage = L"The max method should not throw an exception for two 5x5 matrices"
					". The following was thrown: " + std::wstring(e.what(), e.what() + std::strlen(e.what()));
				Assert::Fail(errorMessage.c_str());

			}
		}


		TEST_METHOD(FloatMatrix_ApplyTest_Valid5x5) {

			// Set up context
			Matrix2D_S<float> fiveXFive_1(5, 5, HORIZONTAL, 1.0f);

			// Assert No Exception
			try {

				Matrix2D_S<float> output = fiveXFive_1.apply([](float a) { return a * 2; });
				float** outputMatrixPtr = output.getMatrix();

				Assert::AreEqual((long) 5, output.getCols(), L"The number of input columns do not match");
				Assert::AreEqual((long) 5, output.getRows(), L"The number of input columns do not match");
				Assert::AreEqual((long) 3, output.getPadding(), L"The amount of padding does not match the correct SIMD quantity");

				for (long i = 0; i < output.getRows(); i++) {
					for (long j = 0; j < output.getCols(); j++) {
						Assert::AreEqual(2, (int) outputMatrixPtr[i][j], L"The function was not applied correctly.");
					}
				}

			} catch (std::invalid_argument& e) {

				std::wstring errorMessage = L"The apply method should not throw an exception for a 5x5 matrice"
					". The following was thrown: " + std::wstring(e.what(), e.what() + std::strlen(e.what()));
				Assert::Fail(errorMessage.c_str());

			}
		}


		TEST_METHOD(FloatMatrixV_ApplyTest_Valid5x5) {

			// Set up context
			Matrix2D_S<float> fiveXFive_1(5, 5, VERTICAL, 1.0f);

			// Assert No Exception
			try {

				Matrix2D_S<float> output = fiveXFive_1.apply([](float a) { return a * 2; });
				float** outputMatrixPtr = output.getMatrix();

				Assert::AreEqual((long) 5, output.getCols(), L"The number of input columns do not match");
				Assert::AreEqual((long) 5, output.getRows(), L"The number of input columns do not match");
				Assert::AreEqual((long) 3, output.getPadding(), L"The amount of padding does not match the correct SIMD quantity");

				for (long i = 0; i < output.getRows(); i++) {
					for (long j = 0; j < output.getCols(); j++) {
						Assert::AreEqual(2, (int) outputMatrixPtr[i][j], L"The function was not applied correctly.");
					}
				}

			} catch (std::invalid_argument& e) {

				std::wstring errorMessage = L"The apply method should not throw an exception for a 5x5 matrice"
					". The following was thrown: " + std::wstring(e.what(), e.what() + std::strlen(e.what()));
				Assert::Fail(errorMessage.c_str());

			}
		}


		TEST_METHOD(FloatMatrix_TransposeTest_Valid5x5) {

			// Set up context
			Matrix2D_S<float> fiveXFive_1(5, 1, HORIZONTAL, 1.0f);

			// Assert No Exception
			try {

				Matrix2D_S<float> output = fiveXFive_1.transpose();

				Assert::AreEqual((long)5, output.getCols(), L"The number of input columns do not match");
				Assert::AreEqual((long)1, output.getRows(), L"The number of input columns do not match");
				Assert::AreEqual((long)3, output.getPadding(), L"The amount of padding does not match the correct SIMD quantity");


			} catch (std::invalid_argument& e) {

				std::wstring errorMessage = L"The transpose method should not throw an exception for a 5x1 matrice"
					". The following was thrown: " + std::wstring(e.what(), e.what() + std::strlen(e.what()));
				Assert::Fail(errorMessage.c_str());

			}
		}


		TEST_METHOD(FloatMatrixV_TransposeTest_Valid5x5) {

			// Set up context
			Matrix2D_S<float> fiveXFive_1(5, 1, VERTICAL, 1.0f);

			// Assert No Exception
			try {

				Matrix2D_S<float> output = fiveXFive_1.transpose();
				float** outputMatrixPtr = output.getMatrix();

				Assert::AreEqual((long) 5, output.getCols(), L"The number of input columns do not match");
				Assert::AreEqual((long) 1, output.getRows(), L"The number of input columns do not match");
				Assert::AreEqual((long) 7, output.getPadding(), L"The amount of padding does not match the correct SIMD quantity");

			} catch (std::invalid_argument& e) {

				std::wstring errorMessage = L"The transpose method should not throw an exception for a 5x1 matrice"
					". The following was thrown: " + std::wstring(e.what(), e.what() + std::strlen(e.what()));
				Assert::Fail(errorMessage.c_str());

			}
		}


		TEST_METHOD(FloatMatrix_RealignTestHtoV_Valid5x1) {

			// Set up context
			Matrix2D_S<float> fiveXFive_2(5, 1, HORIZONTAL, 1.0f);

			// Assert No Exception
			try {

				Matrix2D_S<float> output_1 = fiveXFive_2.realign();
				float** outputMatrixPtr = output_1.getMatrix();

				for (long i = 0; i < output_1.getCols(); i++) {
					for (long j = 0; j < output_1.getRows(); j++) {
						Assert::AreEqual(1, (int)outputMatrixPtr[i][j], L"The re-alignment was not performed correctly");
					}
				}

				Assert::AreEqual((long) 1, output_1.getCols(), L"The number of input columns do not match");
				Assert::AreEqual((long) 5, output_1.getRows(), L"The number of input columns do not match");
				Assert::AreEqual((long) 3, output_1.getPadding(), L"The amount of padding does not match the correct SIMD quantity");

			} catch (std::invalid_argument& e) {

				std::wstring errorMessage = L"The realign method should not throw an exception for a 5x1 matrice"
					". The following was thrown: " + std::wstring(e.what(), e.what() + std::strlen(e.what()));
				Assert::Fail(errorMessage.c_str());

			}
		}



		TEST_METHOD(FloatMatrix_RealignTestVtoH_Valid5x1) {

			// Set up context
			Matrix2D_S<float> fiveXFive_2(5, 1, VERTICAL, 1.0f);

			// Assert No Exception
			try {

				Matrix2D_S<float> output_1 = fiveXFive_2.realign();
				float** outputMatrixPtr = output_1.getMatrix();

				for (long i = 0; i < output_1.getRows(); i++) {
					for (long j = 0; j < output_1.getCols(); j++) {
						Assert::AreEqual(1, (int)outputMatrixPtr[i][j], L"The re-alignment was not performed correctly");
					}
				}

				Assert::AreEqual((long) 1, output_1.getCols(), L"The number of input columns do not match");
				Assert::AreEqual((long) 5, output_1.getRows(), L"The number of input columns do not match");
				Assert::AreEqual((long) 7, output_1.getPadding(), L"The amount of padding does not match the correct SIMD quantity");

			} catch (std::invalid_argument& e) {

				std::wstring errorMessage = L"The realign method should not throw an exception for a 5x1 matrice"
					". The following was thrown: " + std::wstring(e.what(), e.what() + std::strlen(e.what()));
				Assert::Fail(errorMessage.c_str());

			}
		}
	};


	TEST_CLASS(DoubleMatrixTests) {
	public:


		TEST_METHOD(DoubleMatrix_Constructor1Test_ValidVector) {

			try {
				// Set up context
				std::vector<std::vector<double>> matrixData(5, std::vector<double>(5, 0));
				Matrix2D_S<double> fiveXFive(matrixData, HORIZONTAL);

			} catch (std::invalid_argument& e) {

				std::wstring errorMessage = L"The constructor should not throw an exception for a valid 5x5 identity matrix. The following was thrown: "
					+ std::wstring(e.what(), e.what() + std::strlen(e.what()));
				Assert::Fail(errorMessage.c_str());

			}

			// Set up context
			std::vector<std::vector<double>> matrixData(5, std::vector<double>(5, 0));
			matrixData[0][0] = 1;
			matrixData[1][1] = 1;
			matrixData[2][2] = 1;
			matrixData[3][3] = 1;
			matrixData[4][4] = 1;
			Matrix2D_S<double> fiveXFive(matrixData, HORIZONTAL);
			double** matrixPtr = fiveXFive.getMatrix();

			// Check if all values are valid
			Assert::AreEqual((long) 5, fiveXFive.getCols(), L"The number of input columns do not match");
			Assert::AreEqual((long) 5, fiveXFive.getRows(), L"The number of input columns do not match");
			Assert::AreEqual((long) 3, fiveXFive.getPadding(), L"The amount of padding does not match the correct SIMD quantity");
			Assert::AreEqual(matrixData[0][0], 1.0, L"The matrix class has not correctly stored the input data");
			Assert::AreEqual(matrixData[1][1], 1.0, L"The matrix class has not correctly stored the input data");
			Assert::AreEqual(matrixData[2][2], 1.0, L"The matrix class has not correctly stored the input data");
			Assert::AreEqual(matrixData[3][3], 1.0, L"The matrix class has not correctly stored the input data");
			Assert::AreEqual(matrixData[4][4], 1.0, L"The matrix class has not correctly stored the input data");
		}


		TEST_METHOD(DoubleMatrix_Constructor1Test_ValidVector_PaddingUnecessary) {

			try {
				// Set up context
				std::vector<std::vector<double>> matrixData(8, std::vector<double>(8, 0));
				Matrix2D_S<double> fiveXFive(matrixData, HORIZONTAL);

			} catch (std::invalid_argument& e) {

				std::wstring errorMessage = L"The constructor should not throw an exception for a valid 5x5 matrix. The following was thrown: "
					+ std::wstring(e.what(), e.what() + std::strlen(e.what()));
				Assert::Fail(errorMessage.c_str());

			}

			// Set up context
			std::vector<std::vector<double>> matrixData(8, std::vector<double>(8, 0));
			Matrix2D_S<double> fiveXFive(matrixData, HORIZONTAL);

			Assert::AreEqual((long) 8, fiveXFive.getCols(), L"The number of input columns do not match");
			Assert::AreEqual((long) 8, fiveXFive.getRows(), L"The number of input columns do not match");
			Assert::AreEqual((long) 0, fiveXFive.getPadding(), L"The amount of padding does not match the correct SIMD quantity");
		}


		TEST_METHOD(DoubleMatrix_Constructor1Test_InvalidVector_InconsistentRowLength) {

			try {
				// Set up context
				std::vector<std::vector<double>> matrixData = {
					std::vector<double>(8, 0),
					std::vector<double>(8, 0),
					std::vector<double>(8, 0),
					std::vector<double>(8, 0),
					std::vector<double>(8, 0),
					std::vector<double>(8, 0),
					std::vector<double>(7, 0), // Incorrect row length
					std::vector<double>(8, 0)
				};
				Matrix2D_S<double> fiveXFive(matrixData, HORIZONTAL);

			} catch (std::invalid_argument& e) {

				Assert::AreEqual(std::wstring(e.what(), e.what() + std::strlen(e.what())), std::wstring(L"Dimension inconsistencies were found in the input data."));

			}
		}


		TEST_METHOD(DoubleMatrix_Constructor2Test_Valid) {

			try {

				Matrix2D_S<double> fiveXFive(5, 5, HORIZONTAL);

			} catch (std::invalid_argument& e) {

				std::wstring errorMessage = L"The constructor should not throw an exception for a valid 5x5. The following was thrown: "
					+ std::wstring(e.what(), e.what() + std::strlen(e.what()));
				Assert::Fail(errorMessage.c_str());

			}

			// Set up context
			Matrix2D_S<double> fiveXFive(5, 5, HORIZONTAL);

			Assert::AreEqual((long)5, fiveXFive.getCols(), L"The number of input columns do not match");
			Assert::AreEqual((long)5, fiveXFive.getRows(), L"The number of input columns do not match");
			Assert::AreEqual((long)3, fiveXFive.getPadding(), L"The amount of padding does not match the correct SIMD quantity");
		}

		TEST_METHOD(DoubleMatrix_CopyConstructorTest_Valid) {

			try {

				Matrix2D_S<double> fiveXFive(5, 5, HORIZONTAL);
				Matrix2D_S<double> fiveXFiveCopy(fiveXFive);

			} catch (std::invalid_argument& e) {

				std::wstring errorMessage = L"The constructor should not throw an exception for a valid 5x5. The following was thrown: "
					+ std::wstring(e.what(), e.what() + std::strlen(e.what()));
				Assert::Fail(errorMessage.c_str());

			}

			// Set up context
			Matrix2D_S<double> fiveXFive(5, 5, HORIZONTAL);
			double** matrixPtr = fiveXFive.getMatrix();
			matrixPtr[0][0] = 1;
			Matrix2D_S<double> fiveXFiveCopy(fiveXFive);
			double** matrixCopyPtr = fiveXFiveCopy.getMatrix();

			Assert::AreEqual((long)5, fiveXFive.getCols(), L"The number of input columns do not match");
			Assert::AreEqual((long)5, fiveXFive.getRows(), L"The number of input columns do not match");
			Assert::AreEqual((long)3, fiveXFive.getPadding(), L"The amount of padding does not match the correct SIMD quantity");
			Assert::AreEqual(matrixCopyPtr[0][0], 1.0, L"The copy didn't correctly copy all of the data");
		}


		TEST_METHOD(DoubleMatrix_Constructor3Test_Valid) {

			try {

				Matrix2D_S<double> fiveXFive(5, 5, HORIZONTAL, 4.2);

			} catch (std::invalid_argument& e) {

				std::wstring errorMessage = L"The constructor should not throw an exception for a valid 5x5. The following was thrown: "
					+ std::wstring(e.what(), e.what() + std::strlen(e.what()));
				Assert::Fail(errorMessage.c_str());

			}

			// Set up context
			Matrix2D_S<double> fiveXFive(5, 5, HORIZONTAL, 4.2);
			double** matrixPtr = fiveXFive.getMatrix();

			Assert::AreEqual((long)5, fiveXFive.getCols(), L"The number of input columns do not match");
			Assert::AreEqual((long)5, fiveXFive.getRows(), L"The number of input columns do not match");
			Assert::AreEqual((long)3, fiveXFive.getPadding(), L"The amount of padding does not match the correct SIMD quantity");
			Assert::AreEqual(matrixPtr[0][0], 4.2, L"The copy didn't correctly copy all of the data");
		}


		TEST_METHOD(DoubleMatrix_Constructor4Test_Valid) {

			try {

				Matrix2D_S<double> fiveXFive(5, 5, HORIZONTAL, []() { return 4.2; });

			} catch (std::invalid_argument& e) {

				std::wstring errorMessage = L"The constructor should not throw an exception for a valid 5x5. The following was thrown: "
					+ std::wstring(e.what(), e.what() + std::strlen(e.what()));
				Assert::Fail(errorMessage.c_str());

			}

			// Set up context
			Matrix2D_S<double> fiveXFive(5, 5, HORIZONTAL, 4.2);
			double** matrixPtr = fiveXFive.getMatrix();

			Assert::AreEqual((long)5, fiveXFive.getCols(), L"The number of input columns do not match");
			Assert::AreEqual((long)5, fiveXFive.getRows(), L"The number of input columns do not match");
			Assert::AreEqual((long)3, fiveXFive.getPadding(), L"The amount of padding does not match the correct SIMD quantity");
			Assert::AreEqual(matrixPtr[0][0], 4.2, L"The copy didn't correctly copy all of the data");
		}


		TEST_METHOD(DoubleMatrix_DotHTest_InvalidColumnsAndRows) {

			// Set up context
			Matrix2D_S<double> fiveXFive_1(5, 5, HORIZONTAL, 2.1);
			Matrix2D_S<double> fiveXFive_2(4, 5, VERTICAL, 0.0);
			double** matrixPtr = fiveXFive_2.getMatrix();

			// Assert No Exception
			try {

				Matrix2D_S<double> output = fiveXFive_1.dot_H(fiveXFive_2);

			} catch (std::invalid_argument& e) {

				Assert::AreEqual(std::wstring(e.what(), e.what() + std::strlen(e.what())), std::wstring(L"The number of columns in matrix 1 does not equal the"
					" number of rows in matrix 2."));

			}
		}



		TEST_METHOD(DoubleMatrix_DotHTest_InvalidMemoryAlignment) {

			// Set up context
			Matrix2D_S<double> fiveXFive_1(5, 5, HORIZONTAL, 2.1);
			Matrix2D_S<double> fiveXFive_2(5, 5, HORIZONTAL, 0.0);
			double** matrixPtr = fiveXFive_2.getMatrix();

			// Assert No Exception
			try {

				Matrix2D_S<double> output = fiveXFive_1.dot_H(fiveXFive_2);

			} catch (std::invalid_argument& e) {

				Assert::AreEqual(std::wstring(e.what(), e.what() + std::strlen(e.what())), std::wstring(L"The orientation of matrix 1 must be horizontal, "
					L"and the orientation of matrix 2 must be vertical - for simd optimization."));

			}
		}


		TEST_METHOD(DoubleMatrix_DotHTest_Valid5x5dot5x5) {

			// Set up context
			Matrix2D_S<double> fiveXFive_1(5, 5, HORIZONTAL, 2.1);
			Matrix2D_S<double> fiveXFive_2(5, 5, VERTICAL, 0.0);
			double** matrixPtr = fiveXFive_2.getMatrix();

			matrixPtr[0][0] = 2;
			matrixPtr[1][1] = 2;
			matrixPtr[2][2] = 2;
			matrixPtr[3][3] = 2;
			matrixPtr[4][4] = 2;

			// Assert No Exception
			try {

				Matrix2D_S<double> output = fiveXFive_1.dot_H(fiveXFive_2);
				double** outputMatrixPtr = output.getMatrix();

				Assert::AreEqual((long)5, output.getCols(), L"The number of input columns do not match");
				Assert::AreEqual((long)5, output.getRows(), L"The number of input columns do not match");
				Assert::AreEqual((long)3, output.getPadding(), L"The amount of padding does not match the correct SIMD quantity");

				for (long i = 0; i < output.getRows(); i++) {
					for (long j = 0; j < output.getCols(); j++) {
						Assert::AreEqual(4.2, outputMatrixPtr[i][j], L"The dot product was not calculate correctly");
					}
				}

			} catch (std::invalid_argument& e) {

				std::wstring errorMessage = L"The dot_H method should not throw an exception for two correctly aligned 5x5"
					" matrices. The following was thrown: " + std::wstring(e.what(), e.what() + std::strlen(e.what()));
				Assert::Fail(errorMessage.c_str());

			}
		}

		TEST_METHOD(DoubleMatrix_DotHTest_Valid5x5dot5x1) {

			// Set up context
			Matrix2D_S<double> fiveXFive_1(5, 5, HORIZONTAL, 2.1);
			Matrix2D_S<double> fiveXFive_2(5, 1, VERTICAL, 1.0);
			double** matrixPtr = fiveXFive_2.getMatrix();


			// Assert No Exception
			try {

				Matrix2D_S<double> output = fiveXFive_1.dot_H(fiveXFive_2);
				double** outputMatrixPtr = output.getMatrix();

				Assert::AreEqual((long) 1, output.getCols(), L"The number of input columns do not match");
				Assert::AreEqual((long) 5, output.getRows(), L"The number of input columns do not match");
				Assert::AreEqual((long) 3, output.getPadding(), L"The amount of padding does not match the correct SIMD quantity");

				for (long i = 0; i < output.getRows(); i++) {
					for (long j = 0; j < output.getCols(); j++) {
						Assert::AreEqual(10.5, outputMatrixPtr[i][j], L"The dot product was not calculate correctly");
					}
				}

			} catch (std::invalid_argument& e) {

				std::wstring errorMessage = L"The dot_H method should not throw an exception for two correctly aligned 5x5"
					" matrices. The following was thrown: " + std::wstring(e.what(), e.what() + std::strlen(e.what()));
				Assert::Fail(errorMessage.c_str());

			}
		}


		TEST_METHOD(DoubleMatrix_DotVTest_Valid5x5dot5x5) {

			// Set up context
			Matrix2D_S<double> fiveXFive_1(5, 5, HORIZONTAL, 2.1);
			Matrix2D_S<double> fiveXFive_2(5, 5, VERTICAL, 0.0);
			double** matrixPtr = fiveXFive_2.getMatrix();

			matrixPtr[0][0] = 2;
			matrixPtr[1][1] = 2;
			matrixPtr[2][2] = 2;
			matrixPtr[3][3] = 2;
			matrixPtr[4][4] = 2;

			// Assert No Exception
			try {

				Matrix2D_S<double> output = fiveXFive_1.dot_V(fiveXFive_2);
				double** outputMatrixPtr = output.getMatrix();

				Assert::AreEqual((long)5, output.getCols(), L"The number of input columns do not match");
				Assert::AreEqual((long)5, output.getRows(), L"The number of input columns do not match");
				Assert::AreEqual((long)3, output.getPadding(), L"The amount of padding does not match the correct SIMD quantity");

				for (long i = 0; i < output.getRows(); i++) {
					for (long j = 0; j < output.getCols(); j++) {
						Assert::AreEqual(4.2, outputMatrixPtr[j][i], L"The dot product was not calculate correctly");
					}
				}

			} catch (std::invalid_argument& e) {

				std::wstring errorMessage = L"The dot_V method should not throw an exception for two correctly aligned 5x5"
					" matrices. The following was thrown: " + std::wstring(e.what(), e.what() + std::strlen(e.what()));
				Assert::Fail(errorMessage.c_str());

			}
		}



		TEST_METHOD(DoubleMatrix_DotVTest_Valid5x5dot5x1) {

			// Set up context
			Matrix2D_S<double> fiveXFive_1(5, 5, HORIZONTAL, 2.1);
			Matrix2D_S<double> fiveXFive_2(5, 1, VERTICAL, 1.0);
			double** matrixPtr = fiveXFive_2.getMatrix();


			// Assert No Exception
			try {

				Matrix2D_S<double> output = fiveXFive_1.dot_V(fiveXFive_2);
				double** outputMatrixPtr = output.getMatrix();

				Assert::AreEqual((long)1, output.getCols(), L"The number of input columns do not match");
				Assert::AreEqual((long)5, output.getRows(), L"The number of input columns do not match");
				Assert::AreEqual((long)3, output.getPadding(), L"The amount of padding does not match the correct SIMD quantity");

				for (long i = 0; i < output.getRows(); i++) {
					for (long j = 0; j < output.getCols(); j++) {
						Assert::AreEqual(10.5, outputMatrixPtr[j][i], L"The dot product was not calculate correctly");
					}
				}

			} catch (std::invalid_argument& e) {

				std::wstring errorMessage = L"The dot_V method should not throw an exception for two correctly aligned 5x5"
					" matrices. The following was thrown: " + std::wstring(e.what(), e.what() + std::strlen(e.what()));
				Assert::Fail(errorMessage.c_str());

			}
		}


		TEST_METHOD(DoubleMatrix_MultiplyTest_Valid5x5xScalar) {

			// Set up context
			Matrix2D_S<double> fiveXFive_1(5, 5, HORIZONTAL, 2.0);

			// Assert No Exception
			try {

				Matrix2D_S<double> output = fiveXFive_1.multiply(3.0);
				double** outputMatrixPtr = output.getMatrix();

				Assert::AreEqual((long)5, output.getCols(), L"The number of input columns do not match");
				Assert::AreEqual((long)5, output.getRows(), L"The number of input columns do not match");
				Assert::AreEqual((long)3, output.getPadding(), L"The amount of padding does not match the correct SIMD quantity");

				for (long i = 0; i < output.getRows(); i++) {
					for (long j = 0; j < output.getCols(); j++) {
						Assert::AreEqual(6, (int)outputMatrixPtr[i][j], L"The multiplication was not calculate correctly");
					}
				}

			} catch (std::invalid_argument& e) {

				std::wstring errorMessage = L"The multiply method should not throw an exception for a 5x5 matrices multiplied"
					" by a scalar. The following was thrown: " + std::wstring(e.what(), e.what() + std::strlen(e.what()));
				Assert::Fail(errorMessage.c_str());

			}
		}


		TEST_METHOD(DoubleMatrix_MultiplyTest_Valid5x5x5x5) {

			// Set up context
			Matrix2D_S<double> fiveXFive_1(5, 5, HORIZONTAL, 2.0);
			Matrix2D_S<double> fiveXFive_2(5, 5, HORIZONTAL, 2.0);

			// Assert No Exception
			try {

				Matrix2D_S<double> output = fiveXFive_1.multiply(fiveXFive_2);
				double** outputMatrixPtr = output.getMatrix();

				Assert::AreEqual((long)5, output.getCols(), L"The number of input columns do not match");
				Assert::AreEqual((long)5, output.getRows(), L"The number of input columns do not match");
				Assert::AreEqual((long)3, output.getPadding(), L"The amount of padding does not match the correct SIMD quantity");

				for (long i = 0; i < output.getRows(); i++) {
					for (long j = 0; j < output.getCols(); j++) {
						Assert::AreEqual(4, (int)outputMatrixPtr[i][j], L"The multiplication was not calculate correctly");
					}
				}

			} catch (std::invalid_argument& e) {

				std::wstring errorMessage = L"The multiply method should not throw an exception for a 5x5 matrices multiplied"
					" by a matrix. The following was thrown: " + std::wstring(e.what(), e.what() + std::strlen(e.what()));
				Assert::Fail(errorMessage.c_str());

			}
		}


		TEST_METHOD(DoubleMatrix_DivideTest_Valid5x5xScalar) {

			// Set up context
			Matrix2D_S<double> fiveXFive_1(5, 5, HORIZONTAL, 2.0);

			// Assert No Exception
			try {

				Matrix2D_S<double> output = fiveXFive_1.divide(2.0);
				double** outputMatrixPtr = output.getMatrix();

				Assert::AreEqual((long)5, output.getCols(), L"The number of input columns do not match");
				Assert::AreEqual((long)5, output.getRows(), L"The number of input columns do not match");
				Assert::AreEqual((long)3, output.getPadding(), L"The amount of padding does not match the correct SIMD quantity");

				for (long i = 0; i < output.getRows(); i++) {
					for (long j = 0; j < output.getCols(); j++) {
						Assert::AreEqual(1, (int)outputMatrixPtr[i][j], L"The division was not calculate correctly");
					}
				}

			} catch (std::invalid_argument& e) {

				std::wstring errorMessage = L"The division method should not throw an exception for a 5x5 matrices divided"
					" by a scalar. The following was thrown: " + std::wstring(e.what(), e.what() + std::strlen(e.what()));
				Assert::Fail(errorMessage.c_str());

			}
		}


		TEST_METHOD(DoubleMatrix_DivideTest_Valid5x5x5x5) {

			// Set up context
			Matrix2D_S<double> fiveXFive_1(5, 5, HORIZONTAL, 2.0);
			Matrix2D_S<double> fiveXFive_2(5, 5, HORIZONTAL, 2.0);

			// Assert No Exception
			try {

				Matrix2D_S<double> output = fiveXFive_1.divide(fiveXFive_2);
				double** outputMatrixPtr = output.getMatrix();

				Assert::AreEqual((long)5, output.getCols(), L"The number of input columns do not match");
				Assert::AreEqual((long)5, output.getRows(), L"The number of input columns do not match");
				Assert::AreEqual((long)3, output.getPadding(), L"The amount of padding does not match the correct SIMD quantity");

				for (long i = 0; i < output.getRows(); i++) {
					for (long j = 0; j < output.getCols(); j++) {
						Assert::AreEqual(1, (int)outputMatrixPtr[i][j], L"The division was not calculate correctly");
					}
				}

			} catch (std::invalid_argument& e) {

				std::wstring errorMessage = L"The divide method should not throw an exception for a 5x5 matrices divided"
					" by a matrix. The following was thrown: " + std::wstring(e.what(), e.what() + std::strlen(e.what()));
				Assert::Fail(errorMessage.c_str());

			}
		}


		TEST_METHOD(DoubleMatrix_AddTest_Valid5x5xScalar) {

			// Set up context
			Matrix2D_S<double> fiveXFive_1(5, 5, HORIZONTAL, 2.0);

			// Assert No Exception
			try {

				Matrix2D_S<double> output = fiveXFive_1.add(2.0);
				double** outputMatrixPtr = output.getMatrix();

				Assert::AreEqual((long)5, output.getCols(), L"The number of input columns do not match");
				Assert::AreEqual((long)5, output.getRows(), L"The number of input columns do not match");
				Assert::AreEqual((long)3, output.getPadding(), L"The amount of padding does not match the correct SIMD quantity");

				for (long i = 0; i < output.getRows(); i++) {
					for (long j = 0; j < output.getCols(); j++) {
						Assert::AreEqual(4, (int)outputMatrixPtr[i][j], L"The addition was not calculate correctly");
					}
				}

			} catch (std::invalid_argument& e) {

				std::wstring errorMessage = L"The add method should not throw an exception for a 5x5 matrices divided"
					" by a scalar. The following was thrown: " + std::wstring(e.what(), e.what() + std::strlen(e.what()));
				Assert::Fail(errorMessage.c_str());

			}
		}


		TEST_METHOD(DoubleMatrix_AddTest_Valid5x5x5x5) {

			// Set up context
			Matrix2D_S<double> fiveXFive_1(5, 5, HORIZONTAL, 2.0);
			Matrix2D_S<double> fiveXFive_2(5, 5, HORIZONTAL, 2.0);

			// Assert No Exception
			try {

				Matrix2D_S<double> output = fiveXFive_1.add(fiveXFive_2);
				double** outputMatrixPtr = output.getMatrix();

				Assert::AreEqual((long)5, output.getCols(), L"The number of input columns do not match");
				Assert::AreEqual((long)5, output.getRows(), L"The number of input columns do not match");
				Assert::AreEqual((long)3, output.getPadding(), L"The amount of padding does not match the correct SIMD quantity");

				for (long i = 0; i < output.getRows(); i++) {
					for (long j = 0; j < output.getCols(); j++) {
						Assert::AreEqual(4, (int)outputMatrixPtr[i][j], L"The addition was not calculate correctly");
					}
				}

			} catch (std::invalid_argument& e) {

				std::wstring errorMessage = L"The add method should not throw an exception for a 5x5 matrices divided"
					" by a matrix. The following was thrown: " + std::wstring(e.what(), e.what() + std::strlen(e.what()));
				Assert::Fail(errorMessage.c_str());

			}
		}


		TEST_METHOD(DoubleMatrix_SubtractTest_Valid5x5xScalar) {

			// Set up context
			Matrix2D_S<double> fiveXFive_1(5, 5, HORIZONTAL, 4.0);

			// Assert No Exception
			try {

				Matrix2D_S<double> output = fiveXFive_1.subtract(2.0);
				double** outputMatrixPtr = output.getMatrix();

				Assert::AreEqual((long)5, output.getCols(), L"The number of input columns do not match");
				Assert::AreEqual((long)5, output.getRows(), L"The number of input columns do not match");
				Assert::AreEqual((long)3, output.getPadding(), L"The amount of padding does not match the correct SIMD quantity");

				for (long i = 0; i < output.getRows(); i++) {
					for (long j = 0; j < output.getCols(); j++) {
						Assert::AreEqual(2, (int)outputMatrixPtr[i][j], L"The subtraction was not calculate correctly");
					}
				}

			} catch (std::invalid_argument& e) {

				std::wstring errorMessage = L"The subtract method should not throw an exception for a 5x5 matrices subtracted"
					" by a scalar. The following was thrown: " + std::wstring(e.what(), e.what() + std::strlen(e.what()));
				Assert::Fail(errorMessage.c_str());

			}
		}


		TEST_METHOD(DoubleMatrix_SubtractTest_Valid5x5x5x5) {

			// Set up context
			Matrix2D_S<double> fiveXFive_1(5, 5, HORIZONTAL, 4.0);
			Matrix2D_S<double> fiveXFive_2(5, 5, HORIZONTAL, 2.0);

			// Assert No Exception
			try {

				Matrix2D_S<double> output = fiveXFive_1.subtract(fiveXFive_2);
				double** outputMatrixPtr = output.getMatrix();

				Assert::AreEqual((long)5, output.getCols(), L"The number of input columns do not match");
				Assert::AreEqual((long)5, output.getRows(), L"The number of input columns do not match");
				Assert::AreEqual((long)3, output.getPadding(), L"The amount of padding does not match the correct SIMD quantity");

				for (long i = 0; i < output.getRows(); i++) {
					for (long j = 0; j < output.getCols(); j++) {
						Assert::AreEqual(2, (int)outputMatrixPtr[i][j], L"The subtraction was not calculate correctly");
					}
				}

			} catch (std::invalid_argument& e) {

				std::wstring errorMessage = L"The subtract method should not throw an exception for a 5x5 matrices subtracted"
					" by a matrix. The following was thrown: " + std::wstring(e.what(), e.what() + std::strlen(e.what()));
				Assert::Fail(errorMessage.c_str());

			}
		}


		TEST_METHOD(DoubleMatrix_MinTest_Valid5x5xScalar) {

			// Set up context
			Matrix2D_S<double> fiveXFive_1(5, 5, HORIZONTAL, 4.0);

			// Assert No Exception
			try {

				Matrix2D_S<double> output = fiveXFive_1.min(2.0);
				double** outputMatrixPtr = output.getMatrix();

				Assert::AreEqual((long)5, output.getCols(), L"The number of input columns do not match");
				Assert::AreEqual((long)5, output.getRows(), L"The number of input columns do not match");
				Assert::AreEqual((long)3, output.getPadding(), L"The amount of padding does not match the correct SIMD quantity");

				// Validate that the minimum value was found correctly
				for (long i = 0; i < output.getRows(); i++) {
					for (long j = 0; j < output.getCols(); j++) {
						Assert::AreEqual(2, (int)outputMatrixPtr[i][j], L"The minimum was not found correctly");
					}
				}

				for (long i = 0; i < output.getRows(); i++) {
					for (long j = output.getCols(); j < output.getCols() + output.getPadding(); j++) {
						Assert::AreEqual(0, (int)outputMatrixPtr[i][j], L"The minimum function was applied incorrectly to SIMD padding.");
					}
				}

			} catch (std::invalid_argument& e) {

				std::wstring errorMessage = L"The min method should not throw an exception for a 5x5 matrice and"
					" a scalar. The following was thrown: " + std::wstring(e.what(), e.what() + std::strlen(e.what()));
				Assert::Fail(errorMessage.c_str());

			}
		}


		TEST_METHOD(DoubleMatrix_MinTest_Valid5x5x5x5) {

			// Set up context
			Matrix2D_S<double> fiveXFive_1(5, 5, HORIZONTAL, 4.0);
			Matrix2D_S<double> fiveXFive_2(5, 5, HORIZONTAL, 2.0);

			// Assert No Exception
			try {

				Matrix2D_S<double> output = fiveXFive_1.min(fiveXFive_2);
				double** outputMatrixPtr = output.getMatrix();

				Assert::AreEqual((long)5, output.getCols(), L"The number of input columns do not match");
				Assert::AreEqual((long)5, output.getRows(), L"The number of input columns do not match");
				Assert::AreEqual((long)3, output.getPadding(), L"The amount of padding does not match the correct SIMD quantity");

				for (long i = 0; i < output.getRows(); i++) {
					for (long j = 0; j < output.getCols(); j++) {
						Assert::AreEqual(2, (int)outputMatrixPtr[i][j], L"The minimum was not found correctly");
					}
				}

			} catch (std::invalid_argument& e) {

				std::wstring errorMessage = L"The min method should not throw an exception for two 5x5 matrices"
					". The following was thrown: " + std::wstring(e.what(), e.what() + std::strlen(e.what()));
				Assert::Fail(errorMessage.c_str());

			}
		}


		TEST_METHOD(DoubleMatrix_MaxTest_Valid5x5xScalar) {

			// Set up context
			Matrix2D_S<double> fiveXFive_1(5, 5, HORIZONTAL, 0.0);

			// Assert No Exception
			try {

				Matrix2D_S<double> output = fiveXFive_1.max(2.0);
				double** outputMatrixPtr = output.getMatrix();

				Assert::AreEqual((long) 5, output.getCols(), L"The number of input columns do not match");
				Assert::AreEqual((long) 5, output.getRows(), L"The number of input columns do not match");
				Assert::AreEqual((long) 3, output.getPadding(), L"The amount of padding does not match the correct SIMD quantity");

				// Validate that the minimum value was found correctly
				for (long i = 0; i < output.getRows(); i++) {
					for (long j = 0; j < output.getCols(); j++) {
						Assert::AreEqual(2, (int) outputMatrixPtr[i][j], (std::wstring(L"The maximum was not found correctly. Row = ") + std::to_wstring(i) + std::wstring(L" Column = ") + std::to_wstring(j)).c_str());
					}
				}

				for (long i = 0; i < output.getRows(); i++) {
					for (long j = output.getCols(); j < output.getCols() + output.getPadding(); j++) {
						Assert::AreEqual(0, (int) outputMatrixPtr[i][j], L"The maximum function was applied incorrectly to SIMD padding.");
					}
				}

			} catch (std::invalid_argument& e) {

				std::wstring errorMessage = L"The max method should not throw an exception for a 5x5 matrice and"
					" a scalar. The following was thrown: " + std::wstring(e.what(), e.what() + std::strlen(e.what()));
				Assert::Fail(errorMessage.c_str());

			}
		}


		TEST_METHOD(DoubleMatrix_MaxTest_Valid5x5x5x5) {

			// Set up context
			Matrix2D_S<double> fiveXFive_1(5, 5, HORIZONTAL, 4.0);
			Matrix2D_S<double> fiveXFive_2(5, 5, HORIZONTAL, 2.0);

			// Assert No Exception
			try {

				Matrix2D_S<double> output = fiveXFive_1.max(fiveXFive_2);
				double** outputMatrixPtr = output.getMatrix();

				Assert::AreEqual((long)5, output.getCols(), L"The number of input columns do not match");
				Assert::AreEqual((long)5, output.getRows(), L"The number of input columns do not match");
				Assert::AreEqual((long)3, output.getPadding(), L"The amount of padding does not match the correct SIMD quantity");

				for (long i = 0; i < output.getRows(); i++) {
					for (long j = 0; j < output.getCols(); j++) {
						Assert::AreEqual(4, (int)outputMatrixPtr[i][j], L"The maximum was not found correctly");
					}
				}

			} catch (std::invalid_argument& e) {

				std::wstring errorMessage = L"The max method should not throw an exception for two 5x5 matrices"
					". The following was thrown: " + std::wstring(e.what(), e.what() + std::strlen(e.what()));
				Assert::Fail(errorMessage.c_str());

			}
		}

		TEST_METHOD(DoubleMatrixV_MultiplyTest_Valid5x5xScalar) {

			// Set up context
			Matrix2D_S<double> fiveXFive_1(5, 5, VERTICAL, 2.0);

			// Assert No Exception
			try {

				Matrix2D_S<double> output = fiveXFive_1.multiply(3.0);
				double** outputMatrixPtr = output.getMatrix();

				Assert::AreEqual((long)5, output.getCols(), L"The number of input columns do not match");
				Assert::AreEqual((long)5, output.getRows(), L"The number of input columns do not match");
				Assert::AreEqual((long)3, output.getPadding(), L"The amount of padding does not match the correct SIMD quantity");

				for (long i = 0; i < output.getRows(); i++) {
					for (long j = 0; j < output.getCols(); j++) {
						Assert::AreEqual(6, (int)outputMatrixPtr[i][j], L"The multiplication was not calculate correctly");
					}
				}

			} catch (std::invalid_argument& e) {

				std::wstring errorMessage = L"The multiply method should not throw an exception for a 5x5 matrices multiplied"
					" by a scalar. The following was thrown: " + std::wstring(e.what(), e.what() + std::strlen(e.what()));
				Assert::Fail(errorMessage.c_str());

			}
		}


		TEST_METHOD(DoubleMatrixV_MultiplyTest_Valid5x5x5x5) {

			// Set up context
			Matrix2D_S<double> fiveXFive_1(5, 5, VERTICAL, 2.0);
			Matrix2D_S<double> fiveXFive_2(5, 5, VERTICAL, 2.0);

			// Assert No Exception
			try {

				Matrix2D_S<double> output = fiveXFive_1.multiply(fiveXFive_2);
				double** outputMatrixPtr = output.getMatrix();

				Assert::AreEqual((long)5, output.getCols(), L"The number of input columns do not match");
				Assert::AreEqual((long)5, output.getRows(), L"The number of input columns do not match");
				Assert::AreEqual((long)3, output.getPadding(), L"The amount of padding does not match the correct SIMD quantity");

				for (long i = 0; i < output.getRows(); i++) {
					for (long j = 0; j < output.getCols(); j++) {
						Assert::AreEqual(4, (int)outputMatrixPtr[i][j], L"The multiplication was not calculate correctly");
					}
				}

			} catch (std::invalid_argument& e) {

				std::wstring errorMessage = L"The multiply method should not throw an exception for a 5x5 matrices multiplied"
					" by a matrix. The following was thrown: " + std::wstring(e.what(), e.what() + std::strlen(e.what()));
				Assert::Fail(errorMessage.c_str());

			}
		}


		TEST_METHOD(DoubleMatrixV_DivideTest_Valid5x5xScalar) {

			// Set up context
			Matrix2D_S<double> fiveXFive_1(5, 5, VERTICAL, 2.0);

			// Assert No Exception
			try {

				Matrix2D_S<double> output = fiveXFive_1.divide(2.0);
				double** outputMatrixPtr = output.getMatrix();

				Assert::AreEqual((long)5, output.getCols(), L"The number of input columns do not match");
				Assert::AreEqual((long)5, output.getRows(), L"The number of input columns do not match");
				Assert::AreEqual((long)3, output.getPadding(), L"The amount of padding does not match the correct SIMD quantity");

				for (long i = 0; i < output.getRows(); i++) {
					for (long j = 0; j < output.getCols(); j++) {
						Assert::AreEqual(1, (int)outputMatrixPtr[i][j], L"The division was not calculate correctly");
					}
				}

			} catch (std::invalid_argument& e) {

				std::wstring errorMessage = L"The divide method should not throw an exception for a 5x5 matrices divided"
					" by a scalar. The following was thrown: " + std::wstring(e.what(), e.what() + std::strlen(e.what()));
				Assert::Fail(errorMessage.c_str());

			}
		}


		TEST_METHOD(DoubleMatrixV_DivideTest_Valid5x5x5x5) {

			// Set up context
			Matrix2D_S<double> fiveXFive_1(5, 5, VERTICAL, 2.0);
			Matrix2D_S<double> fiveXFive_2(5, 5, VERTICAL, 2.0);

			// Assert No Exception
			try {

				Matrix2D_S<double> output = fiveXFive_1.divide(fiveXFive_2);
				double** outputMatrixPtr = output.getMatrix();

				Assert::AreEqual((long)5, output.getCols(), L"The number of input columns do not match");
				Assert::AreEqual((long)5, output.getRows(), L"The number of input columns do not match");
				Assert::AreEqual((long)3, output.getPadding(), L"The amount of padding does not match the correct SIMD quantity");

				for (long i = 0; i < output.getRows(); i++) {
					for (long j = 0; j < output.getCols(); j++) {
						Assert::AreEqual(1, (int)outputMatrixPtr[i][j], L"The division was not calculate correctly");
					}
				}

			} catch (std::invalid_argument& e) {

				std::wstring errorMessage = L"The divide method should not throw an exception for a 5x5 matrices divided"
					" by a matrix. The following was thrown: " + std::wstring(e.what(), e.what() + std::strlen(e.what()));
				Assert::Fail(errorMessage.c_str());

			}
		}


		TEST_METHOD(DoubleMatrixV_AddTest_Valid5x5xScalar) {

			// Set up context
			Matrix2D_S<double> fiveXFive_1(5, 5, VERTICAL, 2.0);

			// Assert No Exception
			try {

				Matrix2D_S<double> output = fiveXFive_1.add(2.0);
				double** outputMatrixPtr = output.getMatrix();

				Assert::AreEqual((long)5, output.getCols(), L"The number of input columns do not match");
				Assert::AreEqual((long)5, output.getRows(), L"The number of input columns do not match");
				Assert::AreEqual((long)3, output.getPadding(), L"The amount of padding does not match the correct SIMD quantity");

				for (long i = 0; i < output.getRows(); i++) {
					for (long j = 0; j < output.getCols(); j++) {
						Assert::AreEqual(4, (int)outputMatrixPtr[i][j], L"The addition was not calculate correctly");
					}
				}

			} catch (std::invalid_argument& e) {

				std::wstring errorMessage = L"The add method should not throw an exception for a 5x5 matrices divided"
					" by a scalar. The following was thrown: " + std::wstring(e.what(), e.what() + std::strlen(e.what()));
				Assert::Fail(errorMessage.c_str());

			}
		}


		TEST_METHOD(DoubleMatrixV_AddTest_Valid5x5x5x5) {

			// Set up context
			Matrix2D_S<double> fiveXFive_1(5, 5, VERTICAL, 2.0);
			Matrix2D_S<double> fiveXFive_2(5, 5, VERTICAL, 2.0);

			// Assert No Exception
			try {

				Matrix2D_S<double> output = fiveXFive_1.add(fiveXFive_2);
				double** outputMatrixPtr = output.getMatrix();

				Assert::AreEqual((long)5, output.getCols(), L"The number of input columns do not match");
				Assert::AreEqual((long)5, output.getRows(), L"The number of input columns do not match");
				Assert::AreEqual((long)3, output.getPadding(), L"The amount of padding does not match the correct SIMD quantity");

				for (long i = 0; i < output.getRows(); i++) {
					for (long j = 0; j < output.getCols(); j++) {
						Assert::AreEqual(4, (int)outputMatrixPtr[i][j], L"The addition was not calculate correctly");
					}
				}

			} catch (std::invalid_argument& e) {

				std::wstring errorMessage = L"The add method should not throw an exception for a 5x5 matrices divided"
					" by a matrix. The following was thrown: " + std::wstring(e.what(), e.what() + std::strlen(e.what()));
				Assert::Fail(errorMessage.c_str());

			}
		}


		TEST_METHOD(DoubleMatrixV_SubtractTest_Valid5x5xScalar) {

			// Set up context
			Matrix2D_S<double> fiveXFive_1(5, 5, VERTICAL, 4.0);

			// Assert No Exception
			try {

				Matrix2D_S<double> output = fiveXFive_1.subtract(2.0);
				double** outputMatrixPtr = output.getMatrix();

				Assert::AreEqual((long)5, output.getCols(), L"The number of input columns do not match");
				Assert::AreEqual((long)5, output.getRows(), L"The number of input columns do not match");
				Assert::AreEqual((long)3, output.getPadding(), L"The amount of padding does not match the correct SIMD quantity");

				for (long i = 0; i < output.getRows(); i++) {
					for (long j = 0; j < output.getCols(); j++) {
						Assert::AreEqual(2, (int)outputMatrixPtr[i][j], L"The subtraction was not calculate correctly");
					}
				}

			} catch (std::invalid_argument& e) {

				std::wstring errorMessage = L"The subtract method should not throw an exception for a 5x5 matrices subtracted"
					" by a scalar. The following was thrown: " + std::wstring(e.what(), e.what() + std::strlen(e.what()));
				Assert::Fail(errorMessage.c_str());

			}
		}


		TEST_METHOD(DoubleMatrixV_SubtractTest_Valid5x5x5x5) {

			// Set up context
			Matrix2D_S<double> fiveXFive_1(5, 5, VERTICAL, 4.0);
			Matrix2D_S<double> fiveXFive_2(5, 5, VERTICAL, 2.0);

			// Assert No Exception
			try {

				Matrix2D_S<double> output = fiveXFive_1.subtract(fiveXFive_2);
				double** outputMatrixPtr = output.getMatrix();

				Assert::AreEqual((long)5, output.getCols(), L"The number of input columns do not match");
				Assert::AreEqual((long)5, output.getRows(), L"The number of input columns do not match");
				Assert::AreEqual((long)3, output.getPadding(), L"The amount of padding does not match the correct SIMD quantity");

				for (long i = 0; i < output.getRows(); i++) {
					for (long j = 0; j < output.getCols(); j++) {
						Assert::AreEqual(2, (int)outputMatrixPtr[i][j], L"The subtraction was not calculate correctly");
					}
				}

			} catch (std::invalid_argument& e) {

				std::wstring errorMessage = L"The subtract method should not throw an exception for a 5x5 matrices subtracted"
					" by a matrix. The following was thrown: " + std::wstring(e.what(), e.what() + std::strlen(e.what()));
				Assert::Fail(errorMessage.c_str());

			}
		}


		TEST_METHOD(DoubleMatrixV_MinTest_Valid5x5xScalar) {

			// Set up context
			Matrix2D_S<double> fiveXFive_1(5, 5, VERTICAL, 4.0);

			// Assert No Exception
			try {

				Matrix2D_S<double> output = fiveXFive_1.min(2.0);
				double** outputMatrixPtr = output.getMatrix();

				Assert::AreEqual((long) 5, output.getCols(), L"The number of input columns do not match");
				Assert::AreEqual((long) 5, output.getRows(), L"The number of input columns do not match");
				Assert::AreEqual((long) 3, output.getPadding(), L"The amount of padding does not match the correct SIMD quantity");

				for (long i = 0; i < output.getRows(); i++) {
					for (long j = 0; j < output.getCols(); j++) {
						Assert::AreEqual(2, (int)outputMatrixPtr[i][j], L"The minimum was not found correctly");
					}
				}

				for (long i = 0; i < output.getRows(); i++) {
					for (long j = output.getCols(); j < output.getCols() + output.getPadding(); j++) {
						Assert::AreEqual(0, (int)outputMatrixPtr[i][j], L"The minimum function was applied incorrectly to SIMD padding.");
					}
				}

			} catch (std::invalid_argument& e) {

				std::wstring errorMessage = L"The min method should not throw an exception for a 5x5 matrice and"
					" a scalar. The following was thrown: " + std::wstring(e.what(), e.what() + std::strlen(e.what()));
				Assert::Fail(errorMessage.c_str());

			}
		}


		TEST_METHOD(DoubleMatrixV_MinTest_Valid5x5x5x5) {

			// Set up context
			Matrix2D_S<double> fiveXFive_1(5, 5, VERTICAL, 4.0);
			Matrix2D_S<double> fiveXFive_2(5, 5, VERTICAL, 2.0);

			// Assert No Exception
			try {

				Matrix2D_S<double> output = fiveXFive_1.min(fiveXFive_2);
				double** outputMatrixPtr = output.getMatrix();

				Assert::AreEqual((long)5, output.getCols(), L"The number of input columns do not match");
				Assert::AreEqual((long)5, output.getRows(), L"The number of input columns do not match");
				Assert::AreEqual((long)3, output.getPadding(), L"The amount of padding does not match the correct SIMD quantity");

				for (long i = 0; i < output.getRows(); i++) {
					for (long j = 0; j < output.getCols(); j++) {
						Assert::AreEqual(2, (int)outputMatrixPtr[i][j], L"The minimum was not found correctly");
					}
				}

			} catch (std::invalid_argument& e) {

				std::wstring errorMessage = L"The min method should not throw an exception for two 5x5 matrices"
					". The following was thrown: " + std::wstring(e.what(), e.what() + std::strlen(e.what()));
				Assert::Fail(errorMessage.c_str());

			}
		}


		TEST_METHOD(DoubleMatrixV_MaxTest_Valid5x5xScalar) {

			// Set up context
			Matrix2D_S<double> fiveXFive_1(5, 5, VERTICAL, 0.0);

			// Assert No Exception
			try {

				Matrix2D_S<double> output = fiveXFive_1.max(2.0);
				double** outputMatrixPtr = output.getMatrix();

				Assert::AreEqual((long)5, output.getCols(), L"The number of input columns do not match");
				Assert::AreEqual((long)5, output.getRows(), L"The number of input columns do not match");
				Assert::AreEqual((long)3, output.getPadding(), L"The amount of padding does not match the correct SIMD quantity");

				for (long i = 0; i < output.getRows(); i++) {
					for (long j = 0; j < output.getCols(); j++) {
						Assert::AreEqual(2, (int)outputMatrixPtr[i][j], L"The maximum was not found correctly");
					}
				}

				for (long i = 0; i < output.getRows(); i++) {
					for (long j = output.getCols(); j < output.getCols() + output.getPadding(); j++) {
						Assert::AreEqual(0, (int)outputMatrixPtr[i][j], L"The maximum function was applied incorrectly to SIMD padding.");
					}
				}

			} catch (std::invalid_argument& e) {

				std::wstring errorMessage = L"The max method should not throw an exception for a 5x5 matrice and"
					" a scalar. The following was thrown: " + std::wstring(e.what(), e.what() + std::strlen(e.what()));
				Assert::Fail(errorMessage.c_str());

			}
		}


		TEST_METHOD(DoubleMatrixV_MaxTest_Valid5x5x5x5) {

			// Set up context
			Matrix2D_S<double> fiveXFive_1(5, 5, VERTICAL, 4.0);
			Matrix2D_S<double> fiveXFive_2(5, 5, VERTICAL, 2.0);

			// Assert No Exception
			try {

				Matrix2D_S<double> output = fiveXFive_1.max(fiveXFive_2);
				double** outputMatrixPtr = output.getMatrix();

				Assert::AreEqual((long)5, output.getCols(), L"The number of input columns do not match");
				Assert::AreEqual((long)5, output.getRows(), L"The number of input columns do not match");
				Assert::AreEqual((long)3, output.getPadding(), L"The amount of padding does not match the correct SIMD quantity");

				for (long i = 0; i < output.getRows(); i++) {
					for (long j = 0; j < output.getCols(); j++) {
						Assert::AreEqual(4, (int)outputMatrixPtr[i][j], L"The maximum was not found correctly");
					}
				}

			} catch (std::invalid_argument& e) {

				std::wstring errorMessage = L"The max method should not throw an exception for two 5x5 matrices"
					". The following was thrown: " + std::wstring(e.what(), e.what() + std::strlen(e.what()));
				Assert::Fail(errorMessage.c_str());

			}
		}


		TEST_METHOD(DoubleMatrix_ApplyTest_Valid5x5) {

			// Set up context
			Matrix2D_S<double> fiveXFive_1(5, 5, HORIZONTAL, 1.0);

			// Assert No Exception
			try {

				Matrix2D_S<double> output = fiveXFive_1.apply([](double a) { return a * 2; });
				double** outputMatrixPtr = output.getMatrix();

				Assert::AreEqual((long)5, output.getCols(), L"The number of input columns do not match");
				Assert::AreEqual((long)5, output.getRows(), L"The number of input columns do not match");
				Assert::AreEqual((long)3, output.getPadding(), L"The amount of padding does not match the correct SIMD quantity");

				for (long i = 0; i < output.getRows(); i++) {
					for (long j = 0; j < output.getCols(); j++) {
						Assert::AreEqual(2, (int)outputMatrixPtr[i][j], L"The function was not applied correctly.");
					}
				}

			} catch (std::invalid_argument& e) {

				std::wstring errorMessage = L"The apply method should not throw an exception for a 5x5 matrice"
					". The following was thrown: " + std::wstring(e.what(), e.what() + std::strlen(e.what()));
				Assert::Fail(errorMessage.c_str());

			}
		}


		TEST_METHOD(DoubleMatrixV_ApplyTest_Valid5x5) {

			// Set up context
			Matrix2D_S<double> fiveXFive_1(5, 5, VERTICAL, 1.0);

			// Assert No Exception
			try {

				Matrix2D_S<double> output = fiveXFive_1.apply([](double a) { return a * 2; });
				double** outputMatrixPtr = output.getMatrix();

				Assert::AreEqual((long)5, output.getCols(), L"The number of input columns do not match");
				Assert::AreEqual((long)5, output.getRows(), L"The number of input columns do not match");
				Assert::AreEqual((long)3, output.getPadding(), L"The amount of padding does not match the correct SIMD quantity");

				for (long i = 0; i < output.getRows(); i++) {
					for (long j = 0; j < output.getCols(); j++) {
						Assert::AreEqual(2, (int)outputMatrixPtr[i][j], L"The function was not applied correctly.");
					}
				}

			} catch (std::invalid_argument& e) {

				std::wstring errorMessage = L"The apply method should not throw an exception for a 5x5 matrice"
					". The following was thrown: " + std::wstring(e.what(), e.what() + std::strlen(e.what()));
				Assert::Fail(errorMessage.c_str());

			}
		}


		TEST_METHOD(DoubleMatrix_TransposeTest_Valid5x5) {

			// Set up context
			Matrix2D_S<double> fiveXFive_1(5, 1, HORIZONTAL, 1.0);

			// Assert No Exception
			try {

				Matrix2D_S<double> output = fiveXFive_1.transpose();

				Assert::AreEqual((long)5, output.getCols(), L"The number of input columns do not match");
				Assert::AreEqual((long)1, output.getRows(), L"The number of input columns do not match");
				Assert::AreEqual((long)3, output.getPadding(), L"The amount of padding does not match the correct SIMD quantity");


			} catch (std::invalid_argument& e) {

				std::wstring errorMessage = L"The transpose method should not throw an exception for a 5x1 matrice"
					". The following was thrown: " + std::wstring(e.what(), e.what() + std::strlen(e.what()));
				Assert::Fail(errorMessage.c_str());

			}
		}


		TEST_METHOD(DoubleMatrixV_TransposeTest_Valid5x5) {

			// Set up context
			Matrix2D_S<double> fiveXFive_1(5, 1, VERTICAL, 1.0);

			// Assert No Exception
			try {

				Matrix2D_S<double> output = fiveXFive_1.transpose();
				double** outputMatrixPtr = output.getMatrix();

				Assert::AreEqual((long)5, output.getCols(), L"The number of input columns do not match");
				Assert::AreEqual((long)1, output.getRows(), L"The number of input columns do not match");
				Assert::AreEqual((long)3, output.getPadding(), L"The amount of padding does not match the correct SIMD quantity");

			} catch (std::invalid_argument& e) {

				std::wstring errorMessage = L"The transpose method should not throw an exception for a 5x1 matrice"
					". The following was thrown: " + std::wstring(e.what(), e.what() + std::strlen(e.what()));
				Assert::Fail(errorMessage.c_str());

			}
		}


		TEST_METHOD(DoubleMatrix_RealignTestHtoV_Valid5x1) {

			// Set up context
			Matrix2D_S<double> fiveXFive_2(5, 1, HORIZONTAL, 1.0f);

			// Assert No Exception
			try {

				Matrix2D_S<double> output_1 = fiveXFive_2.realign();
				double** outputMatrixPtr = output_1.getMatrix();

				for (long i = 0; i < output_1.getCols(); i++) {
					for (long j = 0; j < output_1.getRows(); j++) {
						Assert::AreEqual(1, (int) outputMatrixPtr[i][j], L"The re-alignment was not performed correctly");
					}
				}

				Assert::AreEqual((long) 1, output_1.getCols(), L"The number of input columns do not match");
				Assert::AreEqual((long) 5, output_1.getRows(), L"The number of input columns do not match");
				Assert::AreEqual((long) 3, output_1.getPadding(), L"The amount of padding does not match the correct SIMD quantity");

			} catch (std::invalid_argument& e) {

				std::wstring errorMessage = L"The realign method should not throw an exception for a 5x1 matrice"
					". The following was thrown: " + std::wstring(e.what(), e.what() + std::strlen(e.what()));
				Assert::Fail(errorMessage.c_str());

			}
		}



		TEST_METHOD(DoubleMatrix_RealignTestVtoH_Valid5x1) {

			// Set up context
			Matrix2D_S<double> fiveXFive_2(5, 1, VERTICAL, 1.0f);

			// Assert No Exception
			try {

				Matrix2D_S<double> output_1 = fiveXFive_2.realign();
				double** outputMatrixPtr = output_1.getMatrix();

				for (long i = 0; i < output_1.getRows(); i++) {
					for (long j = 0; j < output_1.getCols(); j++) {
						Assert::AreEqual(1, (int) outputMatrixPtr[i][j], L"The re-alignment was not performed correctly");
					}
				}

				Assert::AreEqual((long) 1, output_1.getCols(), L"The number of input columns do not match");
				Assert::AreEqual((long) 5, output_1.getRows(), L"The number of input columns do not match");
				Assert::AreEqual((long) 3, output_1.getPadding(), L"The amount of padding does not match the correct SIMD quantity");

			} catch (std::invalid_argument& e) {

				std::wstring errorMessage = L"The realign method should not throw an exception for a 5x1 matrice"
					". The following was thrown: " + std::wstring(e.what(), e.what() + std::strlen(e.what()));
				Assert::Fail(errorMessage.c_str());

			}
		}
	};
};

