#### 2021.09.14

二十世纪十大算法：

1946：John von Neumann, Stan Ulam：Metropolis algorithm，也叫Monte Carlo method

1947：RAND Company, George Dantzig：Simplex求解线性规划

1950：Magnus Hestenes, Eduard Stiefel, and Cornelius Lanczos，国家标准局数值分析研究所（INANBS）：Krylov子空间迭代

1951：Alston Householder of Oak Ridge National Laboratory：Householder约化，形式化的矩阵计算分解

1957：IBM，John Backus，FORTRAN最优编译器

1957-1961：J.G.F. Francis of Ferranti Ltd.：QR算法（矩阵本征值的稳定算法）,$$n\leq 3000$$

1962：Tony Hoare of Elliott Brothers, Ltd.：Quicksort

1965：James Cooley of the IBM T.J. Watson Research Center and John Tukey of Princeton
University and AT&T Bell Laboratories：FFT快速傅里叶

1977：Helaman Ferguson and Rodney Forcade of Brigham Young University：integer relation detection整数关系探测

1987：Leslie Greengard and Vladimir Rokhlin of Yale University：Fast multipole快速多极算法



#### 2021.09.17

问题的敏感性：

输入数据的扰动对于问题解的影响程度的大小：

解对数据波动不敏感称为敏感良态：well-conditioned；

对数据波动异常敏感称为铭感病态：ill-conditioned

条件数：

$$cond=\frac{||问题解的相对变化量||}{||输入数据的波动程度||}$$

函数求值

$$cond=\Big|\frac{[f(\hat{x})-f(x)]/f(x)}{(\hat{x}-x)/x}\Big|\approx|\frac{xf'(x)}{f(x)}|$$

算法稳定性

算法：计算数学从n个已知的数据$x\in R^n$出发，最终获得m个结果$y\in R^m$的一系列有限多次操作集合

好的方法：1.准确的算法；2.计算资源方面快捷有效；

提高算法稳定性：1.避免大数相消；2.避免上溢、下溢；3.避免大数吃小数；4.尽量减少运算次数