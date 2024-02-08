# Matrix-Maths

**Author:** Sten Healey

**Description:** This library is currently a simple matrix class that implements a set of matrix specific mathematical operations (more will be added as needed/requested). Furthermore, the library is implemented in such a way that enables for easy compiler optimization, in particular vectorization.

**Key Features/Limitations** 

Please be aware that many operations require the matrices to have their memory orientation to be configured correctly. Two matrices that have the same memory orientation cannot perform a dot product together.

Similarly elementwise addition and subtraction require the memory orientation of the two matrices to be identical.

These are design limitations that arose because of the benefits that can be gained from SIMD operations and Parallelization which will likely be implemented in a later update. As it stands some compiler will already optimize this code to utilize SIMD.

**Note:** The memory of a given matrice may be realigned using the realign() method.

**Planned Improvements**

...*_Version 1.0:_ This currently implements a simple version of a matrix class, albeit that it is designed to be highly optimized by compilers.

...*_Version 2.0:_ The matrix class will implement methods that explicitly utilize SIMD commands, utilizing the OpenMP compiler commands & library.

...*_Version 3.0:_ The matrix class will implement methods that explicitly utilize SIMD commands in combination with parallelized code utilizing the OpenMP compiler commands & library.

...*_Version 4.0:_ The matrix class will implement methods that are offloaded onto the GPU and parallelized using OpenCL (with SIMD).

...*_Version (proposed) 5.0:_ The matrix class will implement methods that are parallelized heterogeneously using the GPU and CPU (with SIMD).
