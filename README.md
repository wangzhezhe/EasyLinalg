The linear algorithm is the foundation for computer science, beyond what we leant from the text book, I need to implement some simple linear algorithm from scratch servering for different goals.

This repo contains all aspects I know about the liner algorithm from programing's perspective.

StaticMem dir contains several example that create the vector and matrix by static allocation. This is only used for proof of the concepts. All these algorithm are implemented in c. Which is easily to be modified and adapted to other devices which may need particular function signature such as __host__ and __device__. This only serve for specific case (since the assocaited algorithm only works for specific matrix and vector and the code can not be reused by other matrix). The contains here can be used as the example to test file for other more compilcated cases.

DynamicMem dir contains example that use the dynamic memory allocation, the gsl library is a complete and production ready library for this type

StaticMemTemplate shows how to use cpp template to handle the case to support the dynamic vector and matrix size. We just need to implement everything based on template. If the code is designed to run on CPU, we can use the std library, if we want to run same code on GPU or other customized devices, we do not use the std library.

OtherExampels shows some commonly used exmaple based on python scipy, gsl library and the EIGEN library. These libraries is originally designed to run on cpu and it is hard to be adopted to gpu. Since we need to update all function call and adding some key words such as __device__ before the function call.

### TODO

adding eigen decomposition for 3by3 4by4 and 8by8
split out into different header files. defination, symm matrix operation, unsymm ...

### References

GSL online code
https://github.com/ampl/gsl/tree/master

A good tutorial based on c code for linear algebra

https://www.andreinc.net/2021/01/20/writing-your-own-linear-algebra-matrix-library-in-c