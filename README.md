# Matrix-Maths

**Version:** 1.1

**Author:** Sten Healey

### **Description** 
This is a header only matrix library. It implements a matrix class that is mathematically and SIMD optimized (utilizing AVX2 and FMA). With the latest patch in _version 1.1_ it now utilize a more efficient memory access method that virtualizes a 2D array in a 1D space reducing memory access misses and allowing fully contiguous allocation of memory.

### **Key Features/Limitations** 

Please be aware that many operations require the matrices to have their memory orientation configured correctly. Two matrices that have the same memory orientation cannot perform a dot product together. Similarly elementwise addition and subtraction require the memory orientation of the two matrices to be identical. These are design limitations that arose because of the benefits that can be gained from SIMD operations and parallelization which will likely be implemented in a later update.

Benchmarking on complex operations like the dot product has shown that as the matrix grows in size (more than 2000x2000) it can exceed the performance of a naive implementation by 25-30x. There are further mathematical optimizations for divisions by a scalar (the reciprocal is calculated once and then multiplication is performed).

### **Compiling**
1. Please copy the header Matrix2D_S.h into your project

2. If you are compiling in Microsoft Visual Studio with the Microsoft `cl` compiler please add the compiler flag `/arch:AVX2` under `Configuration Properties >> C/C++ >> Command Line`

3. If you are compiling using a GNU style compiler e.g. gcc please add the compiler flags `-mavx2 -mfma`

4. That's all everything should compile correctly!

> [!Note]
> The memory of a given matrice may be realigned using the realign() method.


### **Example Usage**

```cpp
#include <iostream>
#include <cstdio>
#include <stdlib.h>
#include "Matrix2D_S.h" // SIMD Matrix Class

using namespace std;

int main() {

    // Create matrices
    Matrix2D_S<float> fiveXFive_1(5, 5, M_Orientation::VERTICAL, 4.0f); // Vertically oriented matrix
    Matrix2D_S<float> fiveXFive_2(5, 1, M_Orientation::HORIZONTAL, 2.0f); // Horizontally oriented matrix

    // Realign memory
    Matrix2D_S<float> output_1 = fiveXFive_1.realign();
    Matrix2D_S<float> output_2 = fiveXFive_2.realign();

    // Display the input matrices
    std::cout << output_1.toString() << std::endl;
    std::cout << output_2.toString() << std::endl;

    // Perform a dot product
    Matrix2D_S<float> output_3 = output_1.dot_V(output_2);

    // Display the output
    std::cout << output_3.toString() << std::endl;

    return 0;
}
```

### Automated Test Cases
In each directory I have added corresponding automated unit tests. These are implemented for Microsoft Visual Studio Code so they may require some set up in order to run.

### **Planned Patches**
- _Version 1.2:_ I may attempt to implement tiling however this may see diminishing returns given that in order to implement it an additional nested loop is required. However,
there may be benefits when implementing _version 2.0_ which - depending on benchmarking - could see this patch moved into _version 2.0_.

### **Planned Versions - Current Version: 1.0**

- _Version 1.0:_ The matrix class will implement methods that explicitly utilize SIMD commands, utilizing AVX2 and FMA intrinsics.

- _Version 2.0:_ A new matrix class will implement methods that explicitly utilize SIMD commands in combination with multi-processing.

- _Version 3.0:_ A new matrix matrix class will implement methods that are offloaded onto the GPU and parallelized using OpenCL.
