# Matrix-Maths

**Author:** Sten Healey

### **Description** 
This is a header only matrix library. Currently, version 1.0 implements a matrix class that is mathematically and SIMD optimized (utilizing AVX2 and FMA).

### **Key Features/Limitations** 

Please be aware that many operations require the matrices to have their memory orientation configured correctly. Two matrices that have the same memory orientation cannot perform a dot product together. Similarly elementwise addition and subtraction require the memory orientation of the two matrices to be identical. These are design limitations that arose because of the benefits that can be gained from SIMD operations and parallelization which will likely be implemented in a later update.

Thus far testing has revealed a performance improvement of 5x the default behaviour of C++ for complex scenarios like dot products. There are further mathematical optimizations for divisions by a scalar (the reciprocal is calculated once and then multiplication is performed).

### **Compiling**
1. Please copy the header Matrix2D_S.h into your project

2. If you are compiling in Microsoft Visual Studio with the Microsoft `cl` compiler please add the compiler flag `/arch:AVX2` under `Configuration Properties >> C/C++ >> Command Line`

3. If you are compiling using a GNU style compiler e.g. gcc please add the compiler flags `-mavx2 -mfma`

4. That's all everything should compile correctly!

> [!Note]
> The memory of a given matrice may be realigned using the realign() method.


### **Example Usage**

```cpp
#include "Matrix2D_S.h" // The SIMD Matrix Class

int main(){

    Matrix2D_S<float> fiveXFive_1(5, 5, HORIZONTAL, 1.0f); // Horizontal memory alignment
    Matrix2D_S<float> fiveXFive_2(5, 1, VERTICAL, 2.0f); // Vertical memory alignment
    Matrix2D_S<float> output = fiveXFive_1.dot_H(fiveXFive_2); // Dot product

    return 0;
}
```

### Automated Test Cases
In each directory I have added corresponding automated unit tests. These are implemented for Microsoft Visual Studio Code so they may require some set up in order to run.

### **Planned Patches**
- _Version 1.1:_ I need to fix the memory access issue accessing T** isn't cache friendly and the overhead of missing the cache makes it quicker to implement the matrix with T* and accessing it using something like this:

```cpp
T* matrix = aligned_alloc(32, (row_size + padding) * cols);

for(long i=0; i<rows; i++){
    for(long j=0; j<cols; j++){
        matrix[i*(row_size + padding) + j] = 0; // Successfully access memory and when using SIMD steps of 8 or 4 the overhead of the multiplication becomes negligible. 
    }
}
```

I will likely replace getMatrix() with an operator overload - probably something like this: `&T operator[](long a, long b)`. This will help to prevent users from unintentionally accessing padded memory.

### **Planned Versions - Current Version: 1.0**

- _Version 1.0:_ The matrix class will implement methods that explicitly utilize SIMD commands, utilizing AVX2 and FMA intrinsics.

- _Version 2.0:_ A new matrix class will implement methods that explicitly utilize SIMD commands in combination with parallelized code.

- _Version 3.0:_ A new matrix matrix class will implement methods that are offloaded onto the GPU and parallelized using OpenCL.
