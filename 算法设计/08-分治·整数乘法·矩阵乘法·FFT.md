# 08 分治·整数乘法·矩阵乘法·FFT

## 8.1 整数乘法

**问题**：

- **加法**：给定两个 $n$ 位的整数 $a,b$，求 $a+b$ —— $O(n)$ bit operations；
- **乘法**：给定两个 $n$ 位的整数 $a,b$，求 $a\cdot b$ —— 最直观的竖式乘法是 $\Theta(n^2)$ bit operations。

**有没有更好的办法？** 采用分治策略：每个整数分成高位、低位两部分。

记

$$

x=2^{n/2}\cdot x_1+x_0,\qquad y=2^{n/2}\cdot y_1+y_0
$$

则

$$

\begin{aligned}
xy&=2^{n}\cdot x_1y_1+2^{n/2}\cdot(x_1y_0+x_0y_1)+x_0y_0
\end{aligned}
$$

这是 **4 次** $n/2$ 位乘法，对应

$$

T(n)=4T(n/2)+\Theta(n)\ \Rightarrow\ T(n)=\Theta(n^2)
$$

**没有改进。** 改进的目标是：**能不能把 4 次乘法减少？**

### Karatsuba 技巧：4 次乘法 → 3 次

注意到

$$

(x_1+x_0)(y_1+y_0)=x_1y_1+x_1y_0+x_0y_1+x_0y_0
$$

正好是上面四个乘积之和。于是

$$

xy=2^n\cdot x_1y_1+2^{n/2}\cdot\Big(\underbrace{(x_1+x_0)(y_1+y_0)}_{B}-\underbrace{x_1y_1}_{A}-\underbrace{x_0y_0}_{C}\Big)+x_0y_0
$$

**所以可以把 4 个乘法变成 3 个！**

```
Recursive-Multiply(x, y):
    令 x = 2^{n/2} x1 + x0,  y = 2^{n/2} y1 + y0
    计算 x1 + x0 与 y1 + y0
    p     = Recursive-Multiply(x1 + x0, y1 + y0)
    x1y1  = Recursive-Multiply(x1, y1)
    x0y0  = Recursive-Multiply(x0, y0)
    Return  x1y1 · 2^n + (p - x1y1 - x0y0) · 2^{n/2} + x0y0
```

> **定理 5.13 [Karatsuba–Ofman, 1962]**　Recursive-Multiply 算法在两个 $n$ 位因数上的运行时间是 $O(n^{\log_2 3})=O(n^{1.585})$。

**证明**：

$$

T(n)\le T(\lfloor n/2\rfloor)+T(\lceil n/2\rceil)+T(1+\lceil n/2\rceil)+O(n)
$$

$$

\Rightarrow\ T(n)=O(n^{\log_2 3})=O(n^{1.585})
$$

> **补充**：$\log_2 3\approx1.585$ 来自 [[06-分治·递推关系与逆序计数]] 中的定理 5.4（$T(n)\le qT(n/2)+cn$ 时 $T(n)=O(n^{\log_2 q})$），这里 $q=3$。

---

## 8.2 矩阵乘积

**问题**：给定两个 $n\times n$ 矩阵 $A,B$，求 $C=AB$。

$$

c_{ij}=\sum_{k=1}^{n}a_{ik}b_{kj}
$$

如果采用直接计算，那么复杂度为 $\Theta(n^3)$ 次算术运算。**存不存在更有效的办法？**

### 分块递归

先把 $A,B$ 分成 $\frac n2\times\frac n2$ 的小块：

$$

\begin{bmatrix}C_{11}&C_{12}\\C_{21}&C_{22}\end{bmatrix}
=
\begin{bmatrix}A_{11}&A_{12}\\A_{21}&A_{22}\end{bmatrix}
\times
\begin{bmatrix}B_{11}&B_{12}\\B_{21}&B_{22}\end{bmatrix}
$$

$$

\begin{aligned}
C_{11}&=A_{11}B_{11}+A_{12}B_{21}\\
C_{12}&=A_{11}B_{12}+A_{12}B_{22}\\
C_{21}&=A_{21}B_{11}+A_{22}B_{21}\\
C_{22}&=A_{21}B_{12}+A_{22}B_{22}
\end{aligned}
$$

分别递归计算 **8** 块 $\frac n2\times\frac n2$ 矩阵乘积；组合时把结果通过 **4** 次矩阵加法组合产生。

$$

T(n)=8T(n/2)+\Theta(n^2)\ \Rightarrow\ T(n)=\Theta(n^3)
$$

**没有改进。**

### Strassen 算法：8 次乘法 → 7 次

> **问题的关键在于能否减少乘法的次数。**

$$
\begin{bmatrix}C_{11}&C_{12}\\C_{21}&C_{22}\end{bmatrix}
=
\begin{bmatrix}A_{11}&A_{12}\\A_{21}&A_{22}\end{bmatrix}
\times
\begin{bmatrix}B_{11}&B_{12}\\B_{21}&B_{22}\end{bmatrix}
$$

7 个乘积（每个都只需要 1 次矩阵乘法）：

$$

\begin{aligned}
P_1&=A_{11}\times(B_{12}-B_{22})\\
P_2&=(A_{11}+A_{12})\times B_{22}\\
P_3&=(A_{21}+A_{22})\times B_{11}\\
P_4&=A_{22}\times(B_{21}-B_{11})\\
P_5&=(A_{11}+A_{22})\times(B_{11}+B_{22})\\
P_6&=(A_{12}-A_{22})\times(B_{21}+B_{22})\\
P_7&=(A_{11}-A_{21})\times(B_{11}+B_{12})
\end{aligned}
$$

4 个结果块：

$$

\begin{aligned}
C_{11}&=P_5+P_4-P_2+P_6\\
C_{12}&=P_1+P_2\\
C_{21}&=P_3+P_4\\
C_{22}&=P_5+P_1-P_3-P_7
\end{aligned}
$$

$$

\text{7 次乘法},\qquad 18=10+8 \text{ 次加减法}
$$

> **订正**：课件此处印作"通过 10 次加减法运算，生成**14 个** $\frac n2\times\frac n2$ 矩阵"。按课件自己列出的 7 个式子逐项清点，括号内的和 / 差一共是 **10 个**，与紧跟其后的"$18=10+8$ 次加减法"自洽。故此处取 **10 个**中间矩阵。

### Strassen 算法的结构与复杂度

> **矩阵快速乘积运算 (Strassen, 1969)**
>
> 1. 先把 $A,B$ 分成 $\frac n2\times\frac n2$ 的小块；
> 2. **计算**：通过 **10 次加减法**运算，生成 10 个 $\frac n2\times\frac n2$ 矩阵；
> 3. **分治**：递归生成 **7** 个 $\frac n2\times\frac n2$ 的矩阵乘积；
> 4. **组合**：7 个矩阵乘积通过 **8 次矩阵的加减法**生成需要的四个矩阵乘积。

$$

T(n)=7T(n/2)+\Theta(n^2)\ \Rightarrow\ T(n)=\Theta(n^{\log_2 7})=O(n^{2.81})
$$

> **注意**：这是一个**理论上的界**，$n$ 比较大的时候才能体现出 Strassen 算法的优势；比如 $n\sim 2{,}500$ 时，Strassen 算法在苹果机上的运行速度是常规算法的八十多倍。

### 分块矩阵乘法的下界问题

> - 两个 $2\times2$ 分块矩阵乘法能否用 7 次向量乘法实现？
>   **可以** [Strassen, 1969]，$\Theta(n^{\log_2 7})=O(n^{2.81})$。
> - 再问：两个 $2\times2$ 分块矩阵乘法能否用 $2\times2$ 矩阵乘法…、能否用 6 次向量乘法实现？
>   **不可能** [Hopcroft and Kerr 1971]。
> - 两个 $3\times3$ 分块矩阵乘法能否用 21 次向量乘法实现？
>   **可以**，$\Theta(n^{\log_3 21})=O(n^{2.77})$。
> - 两个 $70\times70$ 分块矩阵乘法能否用 143,640 次向量乘法实现？
>   **可以**，$O(\log_{70}143640)=O(n^{2.80})$。

### 复杂度最低记录（截至课件）

- December 1979：$O(n^{2.521813})$；
- January 1980：$O(n^{2.52181})$；
- 1987，Coppersmith–Winograd：$O(n^{2.376})$；
- 2020 年 10 月，Alman、Vassilevska Williams：$O(n^{2.3728596})$。

> **猜想**：对于任意的 $\varepsilon>0$，存在 $O(n^{2+\varepsilon})$ 的优化算法。

> **补充**：指数越接近 2，算法中的常数越大、**实际可用性就越低**——需要 $n$ 达到一定规模（通常是数千甚至数万）才能体现优势。

---

## 8.3 卷积与快速傅里叶变换 (FFT)

### 应用场景

signal processing、speech recognition、data compression、image processing；DVD、JPEG、MP3；Numerical solutions to Poisson's equation。

> *The FFT is one of the truly great computational developments of 20th century. It has changed the face of science and engineering so much that it is not an exaggeration to say that life as we know it would be very different without the FFT.* —— Charles Van Loan

### 多项式的两种表示

**系数表示法**：

$$

A(x)=a_0+a_1x+a_2x^2+\cdots+a_{n-1}x^{n-1},\qquad B(x)=b_0+b_1x+\cdots+b_{n-1}x^{n-1}
$$

| 运算 | 系数表示的代价 | 说明 |
|---|---|---|
| **求值** | $O(n)$ | 采用 **Horner 方法**（秦九韶）：$A(x)=a_0+x(a_1+x(a_2+\cdots+x(a_{n-2}+x a_{n-1})\cdots))$ |
| **加法** | $O(n)$ | $A(x)+B(x)=(a_0+b_0)+(a_1+b_1)x+\cdots+(a_{n-1}+b_{n-1})x^{n-1}$ |
| **乘法（卷积）** | $O(n^2)$ | $A(x)\times B(x)=\displaystyle\sum_{i=0}^{2n-2}c_i x^i$，其中 $c_i=\displaystyle\sum_{j} a_jb_{i-j}$ |

**点值表示法**：

> **代数学基本定理 [Gauss]**　$n$ 次复数多项式存在 $n$ 个复数根。
>
> **推论**　$n-1$ 次多项式 $A(x)$ 能够被 $n$ 个不同的 $x$ 处的函数值唯一确定。

$$

A(x):\ (x_0,y_0),(x_1,y_1),\dots,(x_{n-1},y_{n-1}),\qquad y_j=A(x_j)
$$

| 运算 | 点值表示的代价 | 说明 |
|---|---|---|
| **加法** | $O(n)$ | $(x_0,y_0+z_0),\dots,(x_{n-1},y_{n-1}+z_{n-1})$ |
| **乘法** | $O(n)$ | 但需要 $2n-1$ 个点 |
| **求值** | $O(n^2)$ | 相当于作代换 |

### 两种表示之间的转换

$$

\underbrace{a_0,a_1,\dots,a_{n-1}}_{\text{系数表达}}
\ \xrightarrow{\ \text{求值/evaluation}\ }\
\underbrace{(x_0,y_0),\dots,(x_{n-1},y_{n-1})}_{\text{点值表达}}
$$

**系数表达 → 点值表达**：给定 $a_0+a_1x+\cdots+a_{n-1}x^{n-1}$，求它在 $n$ 个不同点 $x_0,\dots,x_{n-1}$ 处的值。

$$

\begin{bmatrix}y_0\\y_1\\\vdots\\y_{n-1}\end{bmatrix}
=
\begin{bmatrix}
1&x_0&x_0^2&\cdots&x_0^{n-1}\\
1&x_1&x_1^2&\cdots&x_1^{n-1}\\
\vdots&\vdots&\vdots&&\vdots\\
1&x_{n-1}&x_{n-1}^2&\cdots&x_{n-1}^{n-1}
\end{bmatrix}
\begin{bmatrix}a_0\\a_1\\\vdots\\a_{n-1}\end{bmatrix}
$$

$O(n^2)$ for matrix-vector multiply；**Vandermonde 矩阵可逆当且仅当 $x_i$ 互不相同**。

**点值表达 → 系数表达**：涉及 Vandermonde 矩阵求逆，直接做的代价是 $O(n^2)$（甚至 $O(n^3)$）。

> **目标**：希望能够找到一种新的表达方式，能够在这两种表达之间进行**高效的转换**。

### 用分治加速求值

**想法**：把多项式表达分成奇偶部分。

$$

A(x)=a_0+a_1x+a_2x^2+a_3x^3+a_4x^4+a_5x^5+a_6x^6+a_7x^7
$$

$$

A_{even}(x)=a_0+a_2x+a_4x^2+a_6x^3,\qquad A_{odd}(x)=a_1+a_3x+a_5x^2+a_7x^3
$$

$$

A(x)=A_{even}(x^2)+x\,A_{odd}(x^2),\qquad A(-x)=A_{even}(x^2)-x\,A_{odd}(x^2)
$$

**尝试 1**：考虑特殊点值 $x=\pm1$：

$$

A(1)=A_{even}(1)+1\cdot A_{odd}(1),\qquad A(-1)=A_{even}(1)-1\cdot A_{odd}(1)
$$

如果计算 2 个次数不超过 $n$ 的多项式在**一点**（$x=1$）的值，那么可以确定一个次数不超过 $n$ 的多项式在**两点**（$x=\pm1$）的值。

**尝试 2**：考虑特殊点值 $x=\pm1,\pm i$：

$$

\begin{aligned}
A(1)&=A_{even}(1)+1\cdot A_{odd}(1), & A(-1)&=A_{even}(1)-1\cdot A_{odd}(1)\\
A(i)&=A_{even}(-1)+i\cdot A_{odd}(-1), & A(-i)&=A_{even}(-1)-i\cdot A_{odd}(-1)
\end{aligned}
$$

如果计算 2 个次数不超过 $n$ 的多项式在**两点**（$x=\pm1$ 处的 $A_{even},A_{odd}$）的值，那么可以确定一个次数不超过 $n$ 的多项式在 **4 点**（$x=\pm1,\pm i$）的值。

**关键之处**：选择

$$

x=\omega^0,\omega^1,\dots,\omega^{n-1}
$$

即 $n$ 次单位根。

### 单位根

> **定义**　如果复数 $x$ 满足 $x^n=1$，那么 $x$ 被称为 **$n$ 次单位根**。
>
> $n$ 次单位根可以表示为 $\omega^0,\omega^1,\dots,\omega^{n-1}$，其中

$$

\omega=e^{2\pi i/n}
$$

> 如果 $n$ 是偶数，那么 $n/2$ 次单位根可以表示为 $\nu^0,\nu^1,\dots,\nu^{n/2-1}$，其中

$$

\nu=e^{4\pi i/n}=\omega^2
$$

### DFT 的定义

目标：求 $A(x)=a_0+a_1x+\cdots+a_{n-1}x^{n-1}$ 在 $n$ 次单位根 $\omega^0,\omega^1,\dots,\omega^{n-1}$ 处的值。

$$

\begin{bmatrix}y_0\\y_1\\y_2\\y_3\\\vdots\\y_{n-1}\end{bmatrix}
=
\begin{bmatrix}
1&1&1&1&\cdots&1\\
1&\omega^1&\omega^2&\omega^3&\cdots&\omega^{n-1}\\
1&\omega^2&\omega^4&\omega^6&\cdots&\omega^{2(n-1)}\\
1&\omega^3&\omega^6&\omega^9&\cdots&\omega^{3(n-1)}\\
\vdots&\vdots&\vdots&\vdots&&\vdots\\
1&\omega^{n-1}&\omega^{2(n-1)}&\omega^{3(n-1)}&\cdots&\omega^{(n-1)(n-1)}
\end{bmatrix}
\begin{bmatrix}a_0\\a_1\\a_2\\a_3\\\vdots\\a_{n-1}\end{bmatrix}
$$

左边称为**离散傅里叶变换 (Discrete Fourier transform)**，中间的矩阵称为**傅里叶矩阵** $F_n$。

### 分治求 DFT

分治策略：把多项式分成偶数部分和奇数部分：

$$

A_{odd}(x)=a_1+a_3x+a_5x^2+\cdots+a_{n-1}x^{(n-2)/2}
$$

$$

A(x)=A_{even}(x^2)+x\,A_{odd}(x^2)
$$

**分别计算**：计算多项式 $A_{even}(x)$、$A_{odd}(x)$ 在 $\frac n2$ 次单位根 $\nu^0,\nu^1,\dots,\nu^{n/2-1}$ 处的值。

**组合**（对 $0\le k<n/2$）：

$$

A(\omega^k)=A_{even}(\nu^k)+\omega^k A_{odd}(\nu^k)
$$

$$

A(\omega^{k+n/2})=A_{even}(\nu^k)-\omega^k A_{odd}(\nu^k)
$$

依据是

$$

\omega^{k+n/2}=-\omega^k,\qquad \nu^k=(\omega^k)^2=(\omega^{k+n/2})^2
$$

> **定理**　FFT 算法能够在 $O(n\log n)$ 步骤内计算 $n-1$ 次多项式在每个 $n$ 次单位根处的值（$n$ 是 2 的整数次幂）。

**运行时间**：

$$

T(2n)=2T(n)+O(n)\ \Rightarrow\ T(n)=O(n\log n)
$$

> **补充（位反序）**：FFT 递归到底层的顺序是 $\nu$ 的下标按二进制位反转排列——000, 100, 010, 110, 001, 101, 011, 111 对应的就是 **bit-reversed order**（perfect shuffle）。迭代实现时先按位反序重排输入即可。

### 傅里叶变换的逆变换 (IFFT)

**目标**：给定 $n-1$ 次多项式在 $n$ 个点 $\omega^0,\omega^1,\dots,\omega^{n-1}$ 处的点值 $y_0,\dots,y_{n-1}$，给出满足条件的唯一的系数 $a_0,\dots,a_{n-1}$。

$$

\mathbf{a}=F_n^{-1}\mathbf{y}
$$

> **求和引理**　设 $\omega$ 是 $n$ 次单位根的生成元，那么

$$

\sum_{j=0}^{n-1}\omega^{kj}=
\begin{cases}
n, & \text{if } k\equiv 0 \bmod n\\
0, & \text{otherwise}
\end{cases}
$$

> **命题**　$F_n\times G_n=I_n$，其中 $G_n$ 的第 $(k',j)$ 项为 $\frac1n\omega^{-k'j}$。

**证明**：

$$

(F_nG_n)_{k,k'}=\frac1n\sum_{j=0}^{n-1}\omega^{kj}\omega^{-jk'}
=\frac1n\sum_{j=0}^{n-1}\omega^{(k-k')j}=
\begin{cases}1,& k=k'\\ 0,& \text{otherwise}\end{cases}
$$

> **如何计算傅里叶变换逆变换？** 发现表达形式与前面的算法计算的对象 (DFT) **非常类似**；应用同样的算法，同时注意把 $\omega^{-1}=e^{-2\pi i/n}$ 看成是 $n$ 次单位根的生成元。

> **定理**　IFFT 算法能够在 $O(n\log n)$ 步骤内，根据 $n$ 次单位根的值，计算出 $n-1$ 次多项式的系数（$n$ 是 2 的整数次幂）。

$$

\underbrace{a_0,a_1,\dots,a_{n-1}}_{\text{系数表达}}
\ \xleftarrow{\ O(n\log n)\ \text{（IFFT）}\ }\
\underbrace{(y_0,\dots,y_{n-1})}_{\text{点值表达}}
$$

### 用 FFT 做卷积

给定向量 $a=(a_0,\dots,a_{n-1})$、$b=(b_0,\dots,b_{n-1})$，把它们看成多项式

$$

A(x)=a_0+a_1x+\cdots+a_{n-1}x^{n-1},\qquad B(x)=b_0+b_1x+\cdots+b_{n-1}x^{n-1}
$$

对于多项式 $C(x)=A(x)B(x)$，可以知道 $c=(c_0,\dots,c_{2n-2})$ 就是**卷积** $a*b$。

```
系数表达 ──── FFT O(n log n) ────> 点值表达
                                      │
                              点值乘法 O(n)
                                      │
系数表达 <─── inverse FFT O(n log n) ── 点值表达
```

> **定理 5.15**　使用快速傅里叶变换，可以在 $O(n\log n)$ 时间内计算初始向量 $a$ 和 $b$ 的卷积。

**应用：计算两个 $n$ 位整数的乘法**

对

$$

a=a_{n-1}\cdots a_1a_0,\qquad b=b_{n-1}\cdots b_1b_0
$$

生成两个 $n$ 次多项式 $p,q$；注意到 $a=p(10)$，$b=q(10)$，求

$$

r(10)=a\cdot b
$$

采用 FFT：$O(n\log n)$。

> **推广**：对于**所有的**正整数 $n$，都存在复杂度为 $O(n\log n)$ 的快速傅里叶变换算法。MIT 的 Frigo and Johnson 完成了优化的 FFT 变换的 C 库文件，见 http://www.fftw.org。

---

## 8.4 本章索引

- 递推关系求解 → [[06-分治·递推关系与逆序计数]]
- 卷积 → [[11-动态规划·序列比对]] 中的序列相似度
- 矩阵乘法的下界问题与"猜想"的表述方式 → [[16-NP与计算的难解性·归约与NP类]] 中的 P vs NP
- 快速傅里叶变换 → [[05-贪心·Huffman编码与数据压缩]] 中的数据压缩应用
