# 10 Hilbert 空间 · 正交分解

> 系列：[[01-预备·Zorn引理与向量空间]] · [[02-度量空间与拓扑空间]] · [[03-拓扑向量空间]] · [[04-赋范空间]] · [[05-Banach空间]] · [[06-可分Banach空间]] · [[07-Banach空间·完备化与紧性]] · [[08-Hilbert空间·基本性质]] · [[09-Hilbert空间·最佳逼近]] · [[11-Hilbert空间·正交基]] · [[12-连续线性算子·有界性与算子范数]] · [[13-连续线性算子·连续线性泛函与弱拓扑]] · [[14-连续线性算子·对偶空间的例子]] · [[15-连续线性算子·伴随算子]] · [[16-连续线性算子·投影算子]]
>
> 3.2 解决了「最佳逼近元存在唯一」。3.3 把它推到极致：当 $B=M$ 是**闭子空间**时，刻画从不等式**升级成等式** $(x-y)\perp M$，于是整个空间被劈成两块互相垂直的部分 —— $H=M\oplus M^\perp$。这是 Hilbert 空间理论里被引用最多的一条定理，也是第 4 章「投影算子」与 3.4「正交基」的共同地基。
>
> 三件事：
>
> 1. **正交与正交补**（定义 3.3.1、命题 3.3.3、3.3.4）：勾股定理；正交补总是闭子空间；$S^{\perp\perp\perp}=S^\perp$。
> 2. **正交投影**（定理 3.3.5）：$M$ 闭子空间时，「$y$ 是最佳逼近元」$\iff$「误差向量 $x-y$ 垂直于 $M$」，于是最优问题化为线性方程组。
> 3. **正交分解定理**（定理 3.3.7）：$H=M\oplus M^\perp$ 且 $M=M^{\perp\perp}$；对一般子集则只有 $S^{\perp\perp}=\overline{\operatorname{span}}S$（定理 3.3.10）。

## 3.3 正交分解

### 3.3.1 正交与正交补

**定义 3.3.1** 设 $X$ 是内积空间。

**(a)** $x,y\in X$ 称为**正交**的，若 $\langle x,y\rangle=0$，记作 $x\perp y$。

**(b)** 设 $S\subseteq X$。$x\perp S$ 表示对一切 $y\in S$ 有 $x\perp y$（显然 $0\perp X$）。两个集合 $S_1,S_2\subseteq X$ 称为**正交**的，若对一切 $x\in S_1$、$y\in S_2$ 有 $x\perp y$，记作 $S_1\perp S_2$。若 $S\perp S$，称 $S$ 为**正交集**。

**(c)** 对 $x\in X$ 定义

$$
x^\perp=\{y\in X:\langle x,y\rangle=0\}.
$$

对子集 $S\subseteq X$ 定义

$$
S^\perp=\{y\in X:\langle x,y\rangle=0,\ \forall x\in S\},
$$

称为 $S$ 的**正交补**（orthogonal complement）。再记 $S^{\perp\perp}=(S^\perp)^\perp$。

**注 3.3.2**

**(a)** $x^\perp$ 是 $X$ 的闭子空间。

**(b)** $S^\perp=(\overline S)^\perp$，其中 $\overline S$ 是 $S$ 在 $X$ 中的闭包。

**证明（(b)）** 由 $S\subseteq\overline S$ 与命题 3.3.4(c) 得 $(\overline S)^\perp\subseteq S^\perp$。反之设 $y\in S^\perp$，取 $S$ 中序列 $z_n\to z\in\overline S$，由内积的连续性（定理 3.1.9）

$$
\langle z,y\rangle=\lim_{n\to\infty}\langle z_n,y\rangle=0,
$$

故 $y\in(\overline S)^\perp$。$\square$

**要点（正交补为什么自动是闭子空间）** 两条等价的看法：

- 由 (b)，$S^\perp$ 只看得见 $S$ 的**闭包**，所以 $S^\perp=(\overline S)^\perp$ 必闭；
- 更直接地，注意 $\langle x,y\rangle=0\iff\langle y,x\rangle=0$，故 $x^\perp=\{y:\langle y,x\rangle=0\}$ 是连续线性泛函 $y\mapsto\langle y,x\rangle$ 的零点集，闭；而

$$
S^\perp=\bigcap_{x\in S}x^\perp
$$

是闭集之交。

**补充（正交补是线性子空间）** 若 $y_1,y_2\in S^\perp$，则对 $x\in S$ 有 $\langle x,\lambda y_1+\mu y_2\rangle=\overline\lambda\langle x,y_1\rangle+\overline\mu\langle x,y_2\rangle=0$（这里用到本书内积对第二变量共轭线性）。故 $S^\perp$ 是子空间。

### 3.3.2 勾股定理与正交补的基本性质

**命题 3.3.3（勾股定理）** 若 $x\perp y$，则

$$
\lVert x+y\rVert^2=\lVert x\rVert^2+\lVert y\rVert^2.
$$

当底域 $K=\mathbb R$ 时逆命题成立。

**证明。** 由 $\langle x,y\rangle=\langle y,x\rangle=0$，

$$
\lVert x+y\rVert^2=\langle x+y,x+y\rangle=\langle x,x\rangle+\langle x,y\rangle+\langle y,x\rangle+\langle y,y\rangle=\lVert x\rVert^2+\lVert y\rVert^2 .
$$

反之，$\lVert x+y\rVert^2-\lVert x\rVert^2-\lVert y\rVert^2=2\operatorname{Re}\langle x,y\rangle$，故等式成立时 $\operatorname{Re}\langle x,y\rangle=0$；当 $K=\mathbb R$ 时 $\langle x,y\rangle=\operatorname{Re}\langle x,y\rangle=0$，即 $x\perp y$。$\square$

**要点（逆命题在复空间失效）** 复情形下「$\operatorname{Re}\langle x,y\rangle=0$」弱于「$\langle x,y\rangle=0$」。取 $X=\mathbb C$、$x=1$、$y=i$：$\lVert 1+i\rVert^2=2=\lVert 1\rVert^2+\lVert i\rVert^2$，但 $\langle 1,i\rangle=i\ne0$。所以复空间的勾股定理**只是充分条件**，不是刻画 —— 这也是复几何与实几何第一处真正分道的地方。

**补充（勾股定理可以迭代）** 若 $x_1,\dots,x_n$ 两两正交，则

$$
\Bigl\lVert\sum_{k=1}^nx_k\Bigr\rVert^2=\sum_{k=1}^n\lVert x_k\rVert^2 .
$$

3.4 的 Bessel 不等式与 Parseval 恒等式就是这条的极限版本。

**命题 3.3.4** 设 $S,S_1$ 是内积空间 $X$ 的子集。

**(a)** $S^\perp$ 是 $X$ 的闭子空间，且 $S\cap S^\perp=\{0\}$；

**(b)** $S\subseteq S^{\perp\perp}$；

**(c)** $S\subseteq S_1\implies S_1^\perp\subseteq S^\perp$；

**(d)** $S^\perp=S^{\perp\perp\perp}$。

**证明。**

**(a)** 由 $S^\perp=\bigcap_{x\in S}x^\perp$ 与注 3.3.2(a) 得 $S^\perp$ 是闭子空间。若 $x_0\in S\cap S^\perp$，则 $\langle x_0,x_0\rangle=0$，故 $x_0=0$。

**(b)** 设 $x\in S$。对任意 $y\in S^\perp$，由 $x\perp y$ 得 $\langle y,x\rangle=\overline{\langle x,y\rangle}=0$，即 $x\perp S^\perp$，故 $x\in S^{\perp\perp}$。

**(c)** 显然：对 $y\in S_1^\perp$，因 $S\subseteq S_1$，$y$ 尤其正交于 $S$ 中每个元素，故 $y\in S^\perp$。

**(d)** 把 (b) 用在 $S^\perp$ 上得 $S^\perp\subseteq(S^\perp)^{\perp\perp}=S^{\perp\perp\perp}$。反向：由 $S\subseteq S^{\perp\perp}$ 与 (c) 得 $S^{\perp\perp\perp}=(S^{\perp\perp})^\perp\subseteq S^\perp$。故相等。$\square$

### 3.3.3 闭子空间的最佳逼近 = 正交投影（定理 3.3.5）

**定理 3.3.5** 设 $M$ 是 Hilbert 空间 $H$ 的闭子空间，$x\in H$。则 $x$ 从 $M$ 的最佳逼近元就是那个**唯一**的 $y\in M$ 使 $(x-y)\perp M$。

**证明。**

**（一）** 设 $y$ 是 $x$ 从 $M$ 的最佳逼近元。对任意 $z\in M$、$\lVert z\rVert=1$，向量 $z'=y+\langle x-y,z\rangle z$ 属于 $M$。于是

$$
\lVert x-y\rVert^2\le\lVert x-z'\rVert^2=\bigl\langle(x-y)-\langle x-y,z\rangle z,\ (x-y)-\langle x-y,z\rangle z\bigr\rangle=\lVert x-y\rVert^2-\lvert\langle x-y,z\rangle\rvert^2 .
$$

最后一步的展开是三项：$\lVert x-y\rVert^2$、$-2\lvert\langle x-y,z\rangle\rvert^2$、$+\lvert\langle x-y,z\rangle\rvert^2$。故 $\lvert\langle x-y,z\rangle\rvert=0$，即 $\langle x-y,z\rangle=0$。对一切单位向量 $z\in M$ 成立 ⟹ $(x-y)\perp M$。

**（二）** 反之设 $y_1\in M$ 使 $(x-y_1)\perp M$。对任意 $u\in M$ 有 $(x-y_1)\perp(y_1-u)$，由勾股定理（命题 3.3.3）

$$
\lVert x-u\rVert^2=\bigl\lVert(x-y_1)+(y_1-u)\bigr\rVert^2=\lVert x-y_1\rVert^2+\lVert y_1-u\rVert^2\ge\lVert x-y_1\rVert^2 .
$$

故 $\lVert x-y_1\rVert=\inf\{\lVert x-u\rVert:u\in M\}$，即 $y_1$ 是最佳逼近元；由定理 3.2.1 的唯一性，$y_1=y$。$\square$

**注 3.3.6** 也称 $y$ 是 $x$ 沿闭子空间 $M$ 的**投影**。

**补充（定理 3.3.5 是定理 3.2.2 的加细）** 一般闭凸集上只有**不等式** (2)：$\operatorname{Re}\langle x-\tilde x,b-\tilde x\rangle\le0$。当 $B=M$ 是子空间时它升级为**等式**。理由：$z\in M\implies y\pm z\in M$，把 $b=y+z$ 与 $b=y-z$ 分别代入得 $\operatorname{Re}\langle x-y,z\rangle\le0$ 与 $\ge0$，故实部为 $0$；再把 $z$ 换成 $iz$（复情形，$iM=M$），得虚部也为 $0$。所以

$$
(x-y)\perp M\iff\operatorname{Re}\langle x-y,z\rangle=0\quad(\forall z\in M).
$$

**必须区分的情形**：$M$ 是闭**凸**集而非子空间时（例如闭球），定理 3.3.5 不成立 —— 误差向量不垂直于任何东西，只有定理 3.2.2 的不等式可用。

**补充（可计算性：Gram 方程组）** 定理 3.3.5 的价值在于把「求最优」变成「解线性方程组」。设 $M=\operatorname{span}\{e_1,\dots,e_n\}$、$y=\sum_{j=1}^nc_je_j$，则 $(x-y)\perp M$ 等价于

$$
\sum_{j=1}^nc_j\langle e_j,e_i\rangle=\langle x,e_i\rangle,\qquad i=1,\dots,n,
$$

即 **Gram 方程组** $Gc=b$，其中 $G=(\langle e_j,e_i\rangle)_{i,j}$、$b_i=\langle x,e_i\rangle$。$G$ 是 Hermite 正定矩阵（当 $e_1,\dots,e_n$ 线性无关时），故 $c$ 唯一。若 $\{e_i\}$ 是标准正交系，则 $G=I$，直接得

$$
c_i=\langle x,e_i\rangle .
$$

最小二乘法、正交投影、Fourier 系数公式全部是这一条的化身。

### 3.3.4 正交分解定理（定理 3.3.7）

**定理 3.3.7（正交分解定理）** 若 $M$ 是 Hilbert 空间 $H$ 的闭子空间，则

$$
H=M\oplus M^\perp,\qquad M=M^{\perp\perp}.
$$

**证明。** 设 $x\in H$。$M$ 是闭凸集，由定理 3.2.1 存在唯一的 $y\in M$ 使

$$
\lVert x-y\rVert=\inf\{\lVert x-z\rVert:z\in M\}.
$$

由定理 3.3.5，$(x-y)\perp M$，故

$$
x=y+(x-y)\in M+M^\perp .
$$

又 $M\cap M^\perp=\{0\}$（命题 3.3.4(a)），故 $H=M\oplus M^\perp$（这里 $\oplus$ 指 3.1.8 意义下的直和）。

再证 $M=M^{\perp\perp}$。由命题 3.3.4(b) 有 $M\subseteq M^{\perp\perp}$。反之设 $z\in M^{\perp\perp}$，按已证的分解写 $z=z_1+z_2$，$z_1\in M$、$z_2\in M^\perp$。因 $z$ 与 $z_1$ 都在 $M^{\perp\perp}$ 中、$z_2\in M^\perp$，故

$$
z_2=z-z_1\in M^{\perp\perp}\cap M^\perp=\{0\},
$$

于是 $z=z_1\in M$。故 $M=M^{\perp\perp}$。$\square$

**注 3.3.8** 若 $H=M+N$ 且 $M\subseteq N^\perp$，则 $M=N^\perp$。

**证明。** 已有 $M\subseteq N^\perp$，只需证 $N^\perp\subseteq M$。设 $x\in N^\perp$，按 $H=M+N$ 写 $x=y+z$，$y\in M$、$z\in N$。则

$$
0=\langle x,z\rangle=\langle y,z\rangle+\langle z,z\rangle=0+\lVert z\rVert^2,
$$

其中 $\langle y,z\rangle=0$ 用了 $y\in M\subseteq N^\perp$。故 $z=0$，$x=y\in M$。$\square$

**要点（注 3.3.8 里 $M+N$ 自动就是直和）** 由 $M\subseteq N^\perp$ 得 $M\cap N=\{0\}$（$u\in M\cap N$ 时 $u\perp N\ni u$，故 $u=0$），所以 $H=M+N$ 其实已经是直和，且这个直和还是**正交直和**。注 3.3.8 的用途：只要能把 $H$ 分解成 $M+N$ 且其中一个与另一个正交，就能断定它是正交分解。

**补充（正交分解的三条常用推论）** 记 $P=P_M$ 为沿 $M$ 的正交投影。

- **算子性质**：$P$ 幂等（$P^2=P$）、自伴（$\langle Px,y\rangle=\langle x,Py\rangle$）、$\lVert P\rVert=1$（$M\ne\{0\}$ 时）。反过来，Hilbert 空间上「幂等 + 自伴」的有界线性算子恰是某个闭子空间的正交投影 —— 这是第 4 章 4.5 节的主题。
- **勾股分解**：$\lVert x\rVert^2=\lVert Px\rVert^2+\lVert x-Px\rVert^2$，即

$$
\lVert x\rVert^2=\lVert Px\rVert^2+d(x,M)^2 .
$$

- **$M^{\perp\perp}=M$ 是「闭」的判据**：一般子集只有 $S^{\perp\perp}=\overline{\operatorname{span}}S$（定理 3.3.10），等号成立 $\iff$ $S$ 本身是闭线性子空间。

### 3.3.5 一个具体例子（例 3.3.9）

**例 3.3.9** 在 $L^2[0,1]$ 中取

$$
S=\Bigl\{f\in L^2[0,1]:f(t)=0\ \text{a.e.}\ t\in\Bigl[0,\frac12\Bigr]\Bigr\}.
$$

则

$$
S^\perp=\Bigl\{g\in L^2[0,1]:g(t)=0\ \text{a.e.}\ t\in\Bigl[\frac12,1\Bigr]\Bigr\},
$$

$$
S^{\perp\perp}=\Bigl\{f\in L^2[0,1]:f(t)=0\ \text{a.e.}\ t\in\Bigl[0,\frac12\Bigr]\Bigr\}=S .
$$

**要点** 这就是正交分解 $H=M\oplus M^\perp$ 在函数空间里的样子：$S$ 是「只活在右半区间上」的函数全体，$S^\perp$ 是「只活在左半区间上」的函数全体，而

$$
L^2[0,1]=S\oplus S^\perp .
$$

$S$ 本身是闭子空间（$L^2$ 收敛蕴含几乎处处收敛于子列，故「在 $[0,\frac12]$ 上 a.e. 为 $0$」在 $L^2$ 极限下保持），所以 $S^{\perp\perp}=S$ 与定理 3.3.7 一致。

### 3.3.6 $S^{\perp\perp}$ 是 $S$ 的闭线性包（定理 3.3.10）

**定理 3.3.10** 设 $S$ 是 Hilbert 空间 $H$ 的非空子集，则 $S^{\perp\perp}$ 是 $S$ 的**闭线性包**（closed linear span）。

**证明。** 设 $M=\overline{\operatorname{span}}S$。由命题 3.3.4，$S^{\perp\perp}$ 是 $H$ 的闭子空间且包含 $S$，故由 $M$ 的最小性（$M$ 是所有含 $S$ 的闭子空间之交）得 $M\subseteq S^{\perp\perp}$。

反之 $S\subseteq M$，由命题 3.3.4(c) 得 $S^{\perp\perp}\subseteq M^{\perp\perp}$，再由定理 3.3.7（$M$ 闭）得 $M^{\perp\perp}=M$。故 $S^{\perp\perp}\subseteq M$。两向合并即 $S^{\perp\perp}=M$。$\square$

**要点（为什么一般不能写 $S^{\perp\perp}=S$）** 由定理 3.3.10，$S^{\perp\perp}=S$ 当且仅当 $S$ 本身就是**闭线性子空间**。极端例子：取 $L^2[0,1]$ 中的多项式全体 $S$。由 Weierstrass 逼近定理的 $L^2$ 版本，$S$ 稠密于 $L^2[0,1]$，于是

$$
S^\perp=\{0\},\qquad S^{\perp\perp}=\{0\}^\perp=L^2[0,1]\supsetneq S .
$$

所以「取两次正交补」这个运算**不是恒等**，而是「取闭线性包」。一句话：$S\mapsto S^{\perp\perp}$ 是一个**闭包算子**（$S\subseteq S^{\perp\perp}$、$S^{\perp\perp\perp\perp}=S^{\perp\perp}$、单调）。

**补充（对偶形式：把距离写成范数的最大值）** 由定理 3.3.7 可得下面这条「距离—范数对偶」公式（练习 3）：设 $M$ 是闭子空间、$x_0\in H$，则

$$
\min\{\lVert x_0-z\rVert:z\in M\}=\max\{\lvert\langle x_0,y\rangle\rvert:y\in M^\perp,\ \lVert y\rVert=1\}.
$$

这个形状值得记住：**左端是「$H$ 中的一个距离」，右端是「$H$ 中一个范数的上确界」**。它其实是 Hahn–Banach 对偶

$$
d(x_0,M)=\sup\{\lvert\phi(x_0)\rvert:\phi\in H^*,\ \lVert\phi\rVert\le1,\ \phi|_M=0\}
$$

在 Hilbert 空间的具体实现 —— 把连续线性泛函 $\phi$ 用 Riesz 表示定理写成 $\phi(\cdot)=\langle\cdot,y\rangle$，则 $\phi|_M=0\iff y\in M^\perp$。第 5 章会把这条一般化。

### 3.3.7 练习

**练习 1** 对 Hilbert 空间 $H$ 的任意非空子集 $M$，$\operatorname{span}M$ 稠密于 $H$ 当且仅当 $M^\perp=\{0\}$。

**解。** 注意 $M^\perp=(\operatorname{span}M)^\perp$（一边由命题 3.3.4(c)，另一边由内积对第一变量的线性性）。记 $E=\overline{\operatorname{span}}M$。

$(\Longrightarrow)$ 若 $E=H$，则 $M^\perp=E^\perp=H^\perp=\{0\}$。

$(\Longleftarrow)$ 若 $M^\perp=\{0\}$，则由定理 3.3.10，$E=E^{\perp\perp}=(M^\perp)^\perp=\{0\}^\perp=H$。$\square$

**练习 2** 设 $M_1,M_2$ 是 Hilbert 空间 $H$ 的子空间。

**(a)** 证明 $(M_1+M_2)^\perp=M_1^\perp\cap M_2^\perp$。

**解。** 由内积对第一变量的线性性，

$$
y\in(M_1+M_2)^\perp\iff\langle y,x_1+x_2\rangle=0\quad\forall x_i\in M_i
$$

$\iff$（取 $x_2=0$ 得 $y\perp M_1$；取 $x_1=0$ 得 $y\perp M_2$）$\iff y\in M_1^\perp\cap M_2^\perp$。$\square$

**(b)** 若 $M_1,M_2$ 闭，证明 $(M_1\cap M_2)^\perp=M_1^\perp+M_2^\perp$。

**注意 —— 这一条按原题（只假设 $M_1,M_2$ 闭）不成立。** 见下。

**正确的无条件结论** 对任意子空间 $M_1,M_2$：

$$
(M_1\cap M_2)^\perp=\overline{M_1^\perp+M_2^\perp},
$$

即只对**闭包**成立。理由是两头一套：

- $\supseteq$：由 (a) 与命题 3.3.4(c)，$M_1^\perp+M_2^\perp\subseteq(M_1\cap M_2)^\perp$，右端是闭的，故闭包也含在内。
- $\subseteq$：由 (a) 用在 $M_1^\perp,M_2^\perp$ 上得 $(M_1^\perp+M_2^\perp)^\perp=M_1^{\perp\perp}\cap M_2^{\perp\perp}=M_1\cap M_2$（末步用了 $M_i$ 闭）。两边取 $\perp$，得 $(M_1\cap M_2)^\perp=\overline{M_1^\perp+M_2^\perp}$。

所以原题要成立，得再补一个闭性假设 —— 例如 $M_1+M_2$ 闭（这是经典的等价条件：$M_1+M_2$ 闭与 $M_1^\perp+M_2^\perp$ 闭互为充要）；此时等号成立。**反例**取 $H=\ell^2$ 与 $H'=H\oplus H$：

$$
M_1=\Gamma(T)=\{(x,Tx):x\in H\},\qquad M_2=H\oplus\{0\},
$$

其中 $T=\operatorname{diag}(1,\tfrac12,\tfrac13,\dots)$ 是有界自伴算子（$\Gamma(T)$ 是 $T$ 的图像，因 $T$ 有界而是 $H'$ 的闭子空间）。此时 $M_1\cap M_2=\{0\}$（$Tx=0\Rightarrow x=0$），故 $(M_1\cap M_2)^\perp=H'$；而

$$
M_1^\perp=\{(-Ty,y):y\in H\},\qquad M_2^\perp=\{0\}\oplus H,
$$

（$M_1^\perp$ 的算法：$\langle(x,Tx),(a,b)\rangle_{H'}=0$ 对一切 $x$ 成立 $\iff a+T^*b=0\iff a=-Tb$。）于是

$$
M_1^\perp+M_2^\perp=\bigl\{(-Ty,\,y+c):y,c\in H\bigr\}=\{(a,d)\in H\oplus H:\ a\in\operatorname{im}T\}=\operatorname{im}T\oplus H .
$$

第一坐标被限制在 $\operatorname{im}T=\{a:\sum n^2\lvert a_n\rvert^2<\infty\}$ 内 —— 这是 $\ell^2$ 的**稠密真子空间**（$a=(1,1/2,1/3,\dots)$ 就不属于它，因 $na_n\equiv1\notin\ell^2$）。所以 $M_1^\perp+M_2^\perp\subsetneq H'=(M_1\cap M_2)^\perp$，而 $\overline{M_1^\perp+M_2^\perp}=H'$ 与闭包版本一致。

**练习 3** 设 $M$ 是 Hilbert 空间 $H$ 的闭子空间、$x_0\in H$，证明

$$
\min\{\lVert x_0-z\rVert:z\in M\}=\max\{\lvert\langle x_0,y\rangle\rvert:y\in M^\perp,\ \lVert y\rVert=1\}.
$$

**解。** 由正交分解（定理 3.3.7），$x_0=x_1+x_2$，$x_1\in M$、$x_2\in M^\perp$。由定理 3.3.5，$x_1$ 是 $x_0$ 从 $M$ 的最佳逼近元，故左端 $=\lVert x_2\rVert$。

右端：对任意 $y\in M^\perp$、$\lVert y\rVert=1$，

$$
\langle x_0,y\rangle=\langle x_1,y\rangle+\langle x_2,y\rangle=0+\langle x_2,y\rangle,
$$

由 Cauchy–Schwarz $\lvert\langle x_2,y\rangle\rvert\le\lVert x_2\rVert$；取 $y=x_2/\lVert x_2\rVert$（若 $x_2=0$ 则两边都为 $0$）取到等号。故右端 $=\lVert x_2\rVert$，等式成立。$\square$

**练习 4** 对 $L^2[0,1]$，求 $M=\{f:\int_0^1f(t)\,dt=0\}$ 的正交补，并求 $d(f_0,M)$，其中 $f_0(x)=e^x$。

**解。** **（1）** 常数函数 $g\equiv c$ 显然满足 $\int_0^1fg=c\int_0^1f=0$，故 $\operatorname{span}\{\mathbf 1\}\subseteq M^\perp$。反之设 $g\in M^\perp$，令

$$
f=g-\Bigl(\int_0^1g\Bigr)\mathbf 1 ,
$$

（减去的是**常数**，不是标量本身）。则 $\int_0^1f=\int_0^1g-\int_0^1g=0$，故 $f\in M$，由 $g\perp f$：

$$
0=\int_0^1fg=\int_0^1g^2-\Bigl(\int_0^1g\Bigr)^2,
$$

即 $g$ 的方差为 $0$，故 $g$ 几乎处处为常数。于是

$$
M^\perp=\operatorname{span}\{\mathbf 1\}=\{\text{常数函数}\}.
$$

**（2）** 由正交分解 $f_0=P_Mf_0+P_{M^\perp}f_0$，而 $P_{M^\perp}f_0$ 是沿常数的投影，即**均值**：

$$
P_{M^\perp}f_0=\Bigl(\int_0^1e^x\,dx\Bigr)\mathbf 1=(e-1)\mathbf 1 .
$$

最佳逼近元是 $P_Mf_0=f_0-P_{M^\perp}f_0=e^x-(e-1)$，**距离则是法向分量的长度**：

$$
d(f_0,M)=\bigl\lVert f_0-P_Mf_0\bigr\rVert_2=\bigl\lVert P_{M^\perp}f_0\bigr\rVert_2=e-1\approx1.7183 .
$$

这里用到 $[0,1]$ 上 $\lVert\mathbf 1\rVert_2=1$，故常数值 $c$ 的范数就是 $\lvert c\rvert$。

> **要点（最容易错的一步）** 距离是 $\lVert P_{M^\perp}f\rVert$，**不是** $\lVert P_Mf\rVert$。本题里零均值部分 $e^x-(e-1)$ 的范数是 $\lVert e^x-(e-1)\rVert_2=\sqrt{(4e-e^2-3)/2}\approx0.4920$，但那是 $\lVert P_Mf_0\rVert$ —— 「$f_0$ 落在 $M$ 里的那一半有多大」，与「$f_0$ 离 $M$ 有多远」是两回事。
>
> 二者由勾股联系：$\lVert f_0\rVert^2=\lVert P_Mf_0\rVert^2+\lVert P_{M^\perp}f_0\rVert^2=0.2420+2.9525=3.1945=(e^2-1)/2$。
>
> 换成统计语言：$\operatorname{Var}(f_0)=\lVert f_0-\bar f_0\rVert_2^2=\lVert P_Mf_0\rVert^2$，而 $d(f_0,M)^2=\lvert\bar f_0\rvert^2$。**方差不是距离；距离是均值**（在区间长度为 $1$ 的归一化下）。

**交叉验算。** 用练习 3 的对偶公式：$M^\perp$ 中的单位向量只有 $\pm\mathbf 1$，故

$$
d(f_0,M)=\max\bigl\{\lvert\langle f_0,y\rangle\rvert:\ y\in M^\perp,\ \lVert y\rVert=1\bigr\}=\Bigl\lvert\int_0^1e^x\,dx\Bigr\rvert=e-1 ,
$$

与上式吻合。**同一节里两道练习能互相验算时一定要验** —— 这一处正是靠它才发现的：若把 $0.4920$ 当距离，就与练习 3 给出的 $e-1$ 直接冲突。

**练习 5** 证明 $\langle x,y\rangle=0$ 当且仅当下列任一条成立：

(i) 对一切 $\alpha\in K$，$\lVert x+\alpha y\rVert\ge\lVert x\rVert$；

(ii) 对一切 $\alpha\in K$，$\lVert x+\alpha y\rVert=\lVert x-\alpha y\rVert$。

**解。** 先算两个恒等式。由内积对第一变量线性、对第二变量共轭线性，

$$
\lVert x+\alpha y\rVert^2=\lVert x\rVert^2+2\operatorname{Re}\bigl(\overline\alpha\langle x,y\rangle\bigr)+\lvert\alpha\rvert^2\lVert y\rVert^2,
$$

$$
\lVert x+\alpha y\rVert^2-\lVert x-\alpha y\rVert^2=4\operatorname{Re}\bigl(\overline\alpha\langle x,y\rangle\bigr).
$$

$(\Longrightarrow)$ 若 $\langle x,y\rangle=0$，两式右端分别变成 $\lVert x\rVert^2+\lvert\alpha\rvert^2\lVert y\rVert^2\ge\lVert x\rVert^2$ 与 $0$，故 (i)(ii) 都成立。

$(\Longleftarrow)$ **由 (i)**：若 $y\ne0$，取 $\alpha=-\langle x,y\rangle/\lVert y\rVert^2$，则 $\operatorname{Re}(\overline\alpha\langle x,y\rangle)=-\lvert\langle x,y\rangle\rvert^2/\lVert y\rVert^2$，代入第一式得

$$
\lVert x+\alpha y\rVert^2=\lVert x\rVert^2-\frac{\lvert\langle x,y\rangle\rvert^2}{\lVert y\rVert^2}\ge\lVert x\rVert^2,
$$

故 $\lvert\langle x,y\rangle\rvert^2\le0$，即 $\langle x,y\rangle=0$。**由 (ii)**：对一切 $\alpha$ 有 $\operatorname{Re}(\overline\alpha\langle x,y\rangle)=0$，取 $\alpha=\langle x,y\rangle$ 得 $\lvert\langle x,y\rangle\rvert^2=0$。$\square$

**要点** (i) 的几何意思是「$x$ 是 $\mathbb C y$ 到 $x$ 的最佳逼近」—— 就是把练习 3 用在 $M=\operatorname{span}\{y\}$ 上的特例。

**练习 6** 设 $M$ 是 Hilbert 空间 $H$ 的闭子空间，证明 $H/M$ 与 $M^\perp$ 等距同构，同构由 $T([x])=P^\perp x$ 给出，其中 $P^\perp$ 是到 $M^\perp$ 的正交投影。

**解。** **良定义**：若 $[x]=[x']$，即 $x-x'\in M$，则 $P^\perp(x-x')=0$，故 $P^\perp x=P^\perp x'$。**线性**：$P^\perp$ 线性。**满射**：若 $w\in M^\perp$ 则 $P^\perp w=w$。**等距**：由商范数的定义与定理 3.3.5，

$$
\lVert[x]\rVert_{H/M}=\inf_{m\in M}\lVert x-m\rVert=\lVert x-P_Mx\rVert=\lVert P^\perp x\rVert .
$$

故 $T$ 是（线性）等距同构。$\square$

**练习 7**

**(i)** 设 $\langle x,z\rangle=\langle y,z\rangle$ 对某个稠密子集 $D\subseteq H$ 中一切 $z$ 成立，证明 $x=y$。

**解。** 令 $w=x-y$，则 $\langle z,w\rangle=0$ 对 $z\in D$ 成立（取了共轭）。$\{z:\langle z,w\rangle=0\}$ 是闭集（内积连续），含稠密集 $D$，故等于 $H$。取 $z=w$ 得 $\lVert w\rVert^2=0$。$\square$

**(ii)** 设 $f,g\in L^2[a,b]$ 且 $\int_a^bf(t)t^ndt=\int_a^bg(t)t^ndt$ 对一切 $n\in\mathbb N$ 成立，证明 $f=g$。

**解。** 令 $w=f-g$。由线性性 $\int_a^bw(t)p(t)\,dt=0$ 对一切多项式 $p$ 成立。多项式全体在 $L^2[a,b]$ 中稠密（Weierstrass 定理 + $[a,b]$ 有界），故 $\int_a^bw(t)h(t)\,dt=0$ 对一切 $h\in L^2[a,b]$ 成立（用 (i) 或直接由 Cauchy–Schwarz 与稠密性过渡）。取 $h=w$ 得 $\lVert w\rVert_2=0$，即 $f=g$ 几乎处处。$\square$

**要点（两条练习是同一件事）** (i) 是抽象版、(ii) 是具体版：**「在稠密集上与所有元素配对都给出 $0$」$\implies$「这个元素是 $0$」**。这就是为什么 Fourier 系数唯一确定函数、矩 $\{\int x^ntf\}$ 唯一确定测度。

---

## 订正

1. **定理 3.3.10 证明里的交叉引用写错了。** 原文末句「$S^{\perp\perp}\subseteq M^{\perp\perp}=M$ by Remark 3.3.7」引的是**注 3.3.7**，但注只有 3.3.6 与 3.3.8 —— 提供 $M^{\perp\perp}=M$ 的是**定理 3.3.7**（正交分解定理）。本文已改。
2. **定理 3.3.5 的证明有一个小缺口。** 原文只对满足 $\lVert z\rVert=1$ 的 $z\in M$ 做了论证，却断言「$\langle x-y,z\rangle=0$ for all $z\in M$」。对 $z\ne0$ 应先对 $z/\lVert z\rVert$ 用同一论证（$0$ 是平凡的）。本文已补。
3. **命题 3.3.4(b) 的证明中写 $\langle y,x\rangle=\langle x,y\rangle$，复情形下不成立。** 应为 $\langle y,x\rangle=\overline{\langle x,y\rangle}$。数值同为 $0$ 故结论不受影响，但等式本身只在实空间成立。本文已改。
4. **练习 2 后半按原题不成立**（见该题解答）：需要额外假设「$M_1+M_2$ 闭」。无条件成立的版本是 $(M_1\cap M_2)^\perp=\overline{M_1^\perp+M_2^\perp}$。已附反例。
5. **练习 3 的左端写错了变量。** 原题给的是 $x_0\in H$，右端用的也是 $x_0$，但左端写成 $\min\{\lVert x_1-z\rVert:z\in M\}$ —— $x_1$ 在题中未定义。应为 $\min\{\lVert x_0-z\rVert:z\in M\}$。本文已改。
6. **注 3.3.2(b) 的 $S^\perp=(\overline S)^\perp$ 在文字抽取时丢了上划线**，原 PDF 里 `(S)` 看起来像 $S$ 本身，容易误读成 $S^\perp=S^\perp$ 这种空话。按上下文与证明（用到了内积连续性）只能是闭包。

---

## 下一步

3.4 **正交基**：把 3.3 的「正交补」用到底 —— 用 Zorn 引理造出极大正交系（正交基），证明每个 Hilbert 空间都有正交基，并得到 $x=\sum\langle x,e_\lambda\rangle e_\lambda$、$\lVert x\rVert^2=\sum\lvert\langle x,e_\lambda\rangle\rvert^2$（Parseval 恒等式）。这也是「Hilbert 空间一定有基、而可分 Banach 空间未必有 Schauder 基」（[[06-可分Banach空间]]）这一对比的落点。

---

系列：[[01-预备·Zorn引理与向量空间]] · [[02-度量空间与拓扑空间]] · [[03-拓扑向量空间]] · [[04-赋范空间]] · [[05-Banach空间]] · [[06-可分Banach空间]] · [[07-Banach空间·完备化与紧性]] · [[08-Hilbert空间·基本性质]] · [[09-Hilbert空间·最佳逼近]] · [[11-Hilbert空间·正交基]] · [[12-连续线性算子·有界性与算子范数]] · [[13-连续线性算子·连续线性泛函与弱拓扑]] · [[14-连续线性算子·对偶空间的例子]] · [[15-连续线性算子·伴随算子]] · [[16-连续线性算子·投影算子]]
