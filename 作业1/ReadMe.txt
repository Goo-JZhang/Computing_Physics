LinearEqs.py存放解线性方程组所需的函数，主要使用CholeskySolve(A,b)以及GaussianElimination(A,b)，分别为Cholesky和高斯消元两种方法解线性方程组的函数，
其中A应当是一个二维list，A[i]为第i行的行向量，A[i][j]为矩阵的第i行第j个元素的值，输出一个list，为线性方程组的解。

HilbertMatrix.py需要和LinearEqs.py处与同一目录下，是作业第三题d问的源代码，直接运行，在同一目录下输出HilbertSolution.txt存放运行结果。

Homework4a.cpp是作业第四题a问的源代码，用c++书写，直接在c++编译器上运行即可，在cpp文件同目录下会生成Homework4aResult.txt的文档以储存结果。

Homework4c.py是作业第四题c问的源码，函数SumFormResult(N)输出c问推导级数在[-N+1,N-1]^3的立方体内的求和值，homework4c()函数输出N从1到10的计算结果。