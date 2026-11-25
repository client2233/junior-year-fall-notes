# 07 Banach 空间 · 完备化与紧性

> 系列：[[01-预备·Zorn引理与向量空间]] · [[02-度量空间与拓扑空间]] · [[03-拓扑向量空间]] · [[04-赋范空间]] · [[05-Banach空间]] · [[06-可分Banach空间]] · [[08-Hilbert空间·基本性质]] · [[09-Hilbert空间·最佳逼近]] · [[10-Hilbert空间·正交分解]] · [[11-Hilbert空间·正交基]] · [[12-连续线性算子·有界性与算子范数]] · [[13-连续线性算子·连续线性泛函与弱拓扑]] · [[14-连续线性算子·对偶空间的例子]] · [[15-连续线性算子·伴随算子]] · [[16-连续线性算子·投影算子]]
>
> 本节三件事，各占一段：
>
> 1. **完备性有个级数版** —— 不必逐个数 Cauchy 列，只需检查「绝对收敛 $\Rightarrow$ 收敛」（定理 2.4.4）。
> 2. **不完备的空间可以补全**，且补全唯一（定理 2.4.5）—— 所以「不完备」不是一个缺陷，只是一个未完成的阶段。
> 3. **无限维空间天生不局部紧**（Riesz 引理 + 定理 2.4.8），这是无限维与 $\mathbb R^n$ 分道扬镳的第一处，也是后面一切「弱拓扑」「弱\* 拓扑」出场的根本原因。
>
> 最后附一件工具：**Banach 不动点定理**（定理 2.4.9），它是「存在唯一性」类定理的标准机器。

## 2.4 Banach 空间中的完备性与紧性

### 2.4.1 级数与绝对收敛

**定义 2.4.1（级数的收敛）** 设 $(x_n)_{n\in\mathbb N}$ 是赋范空间 $E$ 中的序列，称 $S_n=\sum_{i=1}^nx_i$ 为级数 $\sum_{i=1}^\infty x_i$ 的**第 $n$ 个部分和**。若 $S_n$ 在 $E$ 中收敛到 $S$，则称级数 $\sum_{i=1}^\infty x_i$ **收敛**于和 $S$，记 $S=\sum_{i=1}^\infty x_i$；若 $(S_n)_{n\in\mathbb N}$ 发散，则称该级数**发散**。

**命题 2.4.2（Cauchy 判据）** 在 Banach 空间 $E$ 中级数 $\sum_{i=1}^\infty x_i$ 收敛，当且仅当对任意 $\varepsilon>0$ 存在自然数 $N$ 使

$$
\Bigl\lVert\sum_{i=n}^{n+p}x_i\Bigr\rVert<\varepsilon,\qquad\forall n>N,\ p\in\mathbb N.
$$

> **补充** 这正是「Cauchy 列在 Banach 空间中收敛」翻译成部分和的语言。注意右端是**任意 $p$（任意长度的一段尾巴）**，不是相邻两项之差 —— 这一点在把判据用于具体级数时最常被记错。

**定义 2.4.3（绝对收敛）** 级数 $\sum_{i=1}^\infty x_i$ 称**绝对收敛**的，若 $\sum_{i=1}^\infty\lVert x_i\rVert$ 在 $\mathbb K$（$=\mathbb R$ 或 $\mathbb C$）中收敛。

**定理 2.4.4（完备性的级数刻画）** 赋范空间 $E$ 完备，当且仅当 $E$ 中每个绝对收敛级数都收敛。

**证明**

**（$\Rightarrow$）** 设 $E$ 完备，$\sum_{i=1}^\infty\lVert x_i\rVert<\infty$。对 $n>N$、$p\in\mathbb N$，

$$
\Bigl\lVert\sum_{i=n}^{n+p}x_i\Bigr\rVert\le\sum_{i=n}^{n+p}\lVert x_i\rVert.
$$

因 $\sum\lVert x_i\rVert$ 在 $\mathbb K$ 中收敛，由命题 2.4.2（用在一维情形）右端可以任意小；于是左边的部分和列是 $E$ 中的 Cauchy 列，由完备性收敛。

**（$\Leftarrow$）** 反设 $E$ 具有性质「每个绝对收敛级数都收敛」，设 $(x_n)_{n\in\mathbb N}$ 是 $E$ 中的 Cauchy 列。因 Cauchy 列有界，存在 $C>0$ 使 $\lVert x_n\rVert\le C$（$n=1,2,\dots$）。

取子列 $(x_{n_k})_{k\ge1}$ 使

$$
\lVert x_{n_k}-x_{n_{k-1}}\rVert<\frac1{2^k},\qquad\forall k\ge2.
$$

（可行：对每个 $k$ 在 Cauchy 列的足够靠后的位置取一项即可。）于是级数

$$
x_{n_1}+(x_{n_2}-x_{n_1})+\dots+(x_{n_k}-x_{n_{k-1}})+\dots
$$

**绝对收敛**，因为 $\sum_{k\ge2}2^{-k}<\infty$。由假设它收敛，即部分和收敛到某 $x\in E$。但第 $m$ 个部分和恰是 $x_{n_m}$，故

$$
x=\lim_{m\to\infty}x_{n_m}.
$$

断言 $x=\lim_{n\to\infty}x_n$。事实上对任意 $\varepsilon>0$，由 Cauchy 性存在 $N_1>0$ 使

$$
\lVert x_{n+p}-x_n\rVert<\frac\varepsilon2,\qquad\forall n>N_1,\ p\in\mathbb N;
$$

又存在 $k_0>0$ 使 $n_{k_0}>N_1$ 且

$$
\lVert x_{n_k}-x\rVert<\frac\varepsilon2,\qquad\forall k\ge k_0.
$$

于是对一切 $n>n_{k_0}$，

$$
\lVert x_n-x\rVert\le\lVert x_n-x_{n_{k_0}}\rVert+\lVert x_{n_{k_0}}-x\rVert<\frac\varepsilon2+\frac\varepsilon2=\varepsilon,
$$

即 $x=\lim_{n\to\infty}x_n$，故 $E$ 完备。$\square$

> **要点** 这条定理把「完备性」这个只能靠 Cauchy 列**逐个**验证的性质，换成了一个可以靠**比值判别法/比较判别法**验证的性质。后面证明 $L^p$、$\ell^p$ 完备（[[05-Banach空间]] 例 2.2.7）用的正是这条捷径的变体。

### 2.4.2 完备化

**定理 2.4.5（完备化的存在唯一性）** 设 $E$ 是赋范空间。则总存在 Banach 空间 $\widetilde E$ 与从 $E$ 到 $\widetilde E$ 的**等距** $T$，使 $T E$ 在 $\widetilde E$ 中稠密 —— 换言之 $\widetilde E$ 是 $E$ 的**完备化**。并且 $E$ 的完备化在**等距同构**意义下唯一。

**证明（概要）**

在 $E$ 的一切 Cauchy 列构成的向量空间 $\mathcal E$ 上定义等价关系

$$
(x_n)_{n\in\mathbb N}\sim(y_n)_{n\in\mathbb N}\iff\lim_{n\to\infty}\lVert x_n-y_n\rVert=0.
$$

在商空间 $\widetilde E:=\mathcal E/\!\sim$ 上定义

$$
\lVert[(x_n)_{n\in\mathbb N}]\rVert:=\lim_{n\to\infty}\lVert x_n\rVert .
$$

（这个极限存在：$\bigl|\lVert x_n\rVert-\lVert x_m\rVert\bigr|\le\lVert x_n-x_m\rVert$，故 $(\lVert x_n\rVert)$ 是 $\mathbb R$ 中的 Cauchy 列；且它与等价类中代表元的选取无关。）可逐条验证 $\widetilde E$ 在此商范数下是 Banach 空间。

另一方面，映射

$$
x\mapsto[(x,x,x,\dots)]
$$

是 $E$ 到 $\widetilde E$ 的线性等距，其像在 $\widetilde E$ 中稠密。$\square$

> **补充（唯一性怎么证）** 设 $\widetilde E_1$、$\widetilde E_2$ 都是 $E$ 的完备化，$T_i:E\to\widetilde E_i$ 为等距、像稠密。在 $T_1E$ 上定义 $S(T_1x):=T_2x$，则 $S$ 是 $T_1E\to\widetilde E_2$ 的等距。$\widetilde E_2$ 完备，故 $S$ 按**一致连续性**（等距当然一致连续）唯一延拓到 $\overline{T_1E}=\widetilde E_1$ 上，得到一个等距满射 $\widetilde E_1\to\widetilde E_2$。对称地可作逆，故二者等距同构。$\square$
>
> 这个「先在有稠密定义域的等距上定义，再用完备性延拓」的论证，是泛函分析里最常用的三板斧之一。

**例 2.4.6**

(a) 给 $C[a,b]$ 配上 $p$-范数 $\lVert\cdot\rVert_p$（$1\le p<\infty$），则 $(C[a,b],\lVert\cdot\rVert_p)$ **不完备**，其完备化是 $L^p[a,b]$。

(b) $c_{00}$ 的完备化是 $c_0$。

(c) $\mathcal P[0,1]$（多项式空间）的完备化是 $C[0,1]$。

> **补充（怎么读这三条）** 三例都是「**同一批元素、换个范数就换个空间**」的范例：
>
> - (a) 是 [[05-Banach空间]] 例 2.2.2 的直接推广：$\lVert\cdot\rVert_p$ 看不见「单点跳变」，补全时会把这些跳变连成的阶梯函数收进来，收完正好是 $L^p[a,b]$；
> - (b) $\lVert\cdot\rVert_\infty$ 下 $c_{00}$ 的极限恰是「尾巴趋于 $0$」的列，这就是 $c_0$（见 [[05-Banach空间]] 例 2.2.12）；
> - (c) 由 Weierstrass 定理，多项式在 $C[0,1]$ 中稠密而本身不完备，故其完备化就是 $C[0,1]$。
>
> 一句话：**完备化 = 把「按这个范数该收而未收的极限点」补进去。**

### 2.4.3 紧性：Riesz 引理

有限维赋范空间中，闭且有界 $\Rightarrow$ 紧（Bolzano–Weierstrass；结合 [[04-赋范空间]] 定理 2.1.19：$n$ 维赋范空间与 $\mathbb K^n$ 线性且度量同构）。**无限维时这是错的**，而要看清这一点，先需要一个工具。

**定理 2.4.7（Riesz 引理）** 设 $Y$ 是赋范空间 $E$ 的**闭**子空间且 $Y\ne E$，$0<\alpha<1$。则存在 $x_\alpha\in E$，$\lVert x_\alpha\rVert=1$，使

$$
\lVert x_\alpha-y\rVert>\alpha,\qquad\forall y\in Y.
$$

**证明** 因 $Y\ne E$，取 $x\in E\setminus Y$。因 $Y$ 闭，可证

$$
d:=\inf\{\lVert x-y\rVert:\ y\in Y\}>0.
$$

（若 $d=0$，则存在 $y_n\in Y$ 使 $y_n\to x$，$Y$ 闭给出 $x\in Y$，矛盾。）

由 $d<d\alpha^{-1}$ 及下确界的定义，存在 $y_0\in Y$ 使

$$
\lVert x-y_0\rVert<d\alpha^{-1}.
$$

令

$$
x_\alpha=\frac{x-y_0}{\lVert x-y_0\rVert}.
$$

则 $\lVert x_\alpha\rVert=1$。对任意 $y\in Y$，注意 $y_0+\lVert x-y_0\rVert\,y\in Y$（$Y$ 是子空间），故 $\lVert x-(y_0+\lVert x-y_0\rVert\,y)\rVert\ge d$，于是

$$
\lVert x_\alpha-y\rVert=\frac{1}{\lVert x-y_0\rVert}\bigl\lVert x-(y_0+\lVert x-y_0\rVert\,y)\bigr\rVert>(\alpha d^{-1})\,d=\alpha.
$$

$\square$

> **要点**
>
> 1. **为什么要 $\alpha<1$ 而不能取 $\alpha=1$。** 引理只断言可以「几乎垂直」地选到单位向量，不能保证恰好 $d(x_\alpha,Y)=1$。反例：取 $E=C[0,1]$、$Y=\{f:f(0)=0\}$（闭子空间），则常函数 $1$ 到 $Y$ 的距离恰为 $1$ 且**在 $0$ 处取到**（$f\equiv0$ 就是最近点）；但把 $\alpha$ 换成 $1$，一般空间里可能根本没有取到 $1$ 的元素。
> 2. **$\alpha$ 的取值无关紧要**，重要的是「$\alpha$ 与维数无关」这一**一致性**。定理 2.4.8 只用 $\alpha=3/4$。
> 3. 引理可改述为：**闭真子空间外总存在「几乎垂直于 $Y$」的单位向量。** 这是无限维几何与有限维几何的第一个分水岭。

**定理 2.4.8** 设 $E$ 是**无限维**赋范空间。则单位球 $U_E=\{x\in E:\lVert x\rVert\le1\}$ 与单位球面 $S_E=\{x\in E:\lVert x\rVert=1\}$ **都不紧**。

**证明** 取 $x_1\in S_E$。因 $E$ 无限维，$\operatorname{span}\{x_1\}\subsetneq E$；由 [[05-Banach空间]] 命题 2.2.10，$\operatorname{span}\{x_1\}$ 闭。由定理 2.4.7（取 $\alpha=3/4$）存在 $x_2\in S_E$ 使

$$
\lVert x_2-\alpha x_1\rVert\ge\frac34,\qquad\forall\alpha\in\mathbb K.
$$

同理 $\operatorname{span}\{x_1,x_2\}$ 是 $E$ 的真闭子空间，故存在 $x_3\in S_E$ 使

$$
\lVert x_3-\alpha x_1-\beta x_2\rVert\ge\frac34,\qquad\forall\alpha,\beta\in\mathbb K.
$$

如此继续，得到 $S_E$ 中的序列 $(x_n)_{n\in\mathbb N}$ 使

$$
\lVert x_n-x_m\rVert\ge\frac34\qquad(n\ne m).
$$

（对 $m<n$，把 $x_m$ 看作 $\operatorname{span}\{x_1,\dots,x_{n-1}\}$ 中的元素。）这个序列的**任何**子列都不是 Cauchy 列，故没有收敛子列。因此 $S_E$ 不紧；又 $S_E\subseteq U_E$ 是闭子集，若 $U_E$ 紧则 $S_E$ 紧，故 $U_E$ 也不紧。$\square$

> **补充（这条定理的份量）**
>
> - 等价说法：**无限维赋范空间不是局部紧的**（每个点的每个邻域都不紧）。
> - 紧性在度量空间里等价于序列紧，所以上面的「造出一列两两距离 $\ge3/4$」就足以否证紧性。
> - 这是「弱拓扑」出场的根本动机之一：范数拓扑下单位球永不紧，换成 **weak / weak\*** 拓扑后才可能紧（Banach–Alaoglu 定理），而紧性正是极值存在、极限点存在的保障。这条线索通往第 8 章。
> - 注意 $3/4$ 只是个方便的数：任意 $\alpha\in(0,1)$ 都行。

### 2.4.4 Banach 不动点定理

**定理 2.4.9（Banach 不动点定理）** 设 $K$ 是 Banach 空间 $E$ 的**非空闭**子集，$\varphi:K\to K$ 是**压缩**映射，即存在常数 $0<\rho<1$ 使

$$
\lVert\varphi(x)-\varphi(y)\rVert\le\rho\lVert x-y\rVert,\qquad\forall x,y\in K. \tag{2.12}
$$

则 $\varphi$ 在 $K$ 中有唯一不动点，即存在唯一的 $z_0\in K$ 使 $\varphi(z_0)=z_0$。

**证明** 任取 $x_0\in K$，令 $x_{n+1}=\varphi(x_n)$（$n=0,1,2,\dots$）。则

$$
\lVert x_{n+1}-x_n\rVert=\lVert\varphi(x_n)-\varphi(x_{n-1})\rVert\le\rho\lVert x_n-x_{n-1}\rVert\le\rho^2\lVert x_{n-1}-x_{n-2}\rVert\le\dots\le\rho^n\lVert x_1-x_0\rVert.
$$

若 $m>n$，则

$$
\lVert x_m-x_n\rVert\le\sum_{i=n}^{m-1}\lVert x_{i+1}-x_i\rVert\le\sum_{i=n}^{m-1}\rho^i\lVert x_1-x_0\rVert\le\frac{\rho^n}{1-\rho}\lVert x_1-x_0\rVert.
$$

故 $(x_n)_{n\in\mathbb N}$ 是 $E$ 中的 Cauchy 列。由 $E$ 完备，$x_n\to z_0$ 对某 $z_0\in E$；由 $K$ 闭，$z_0\in K$。不等式 (2.12) 说明 $\varphi$ 连续，故 $\varphi(x_n)\to\varphi(z_0)$。在 $x_{n+1}=\varphi(x_n)$ 两端取极限，得 $z_0=\varphi(z_0)$。

**唯一性。** 设 $u\in K$ 且 $\varphi(u)=u$，则由

$$
\lVert u-z_0\rVert=\lVert\varphi(u)-\varphi(z_0)\rVert\le\rho\lVert u-z_0\rVert
$$

与 $0<\rho<1$ 得 $u=z_0$。$\square$

> **要点（三个条件都不可省）**
>
> | 条件 | 省掉会怎样 |
> |---|---|
> | $E$ 完备 | 迭代列只是 Cauchy 列，可能不收敛 |
> | $K$ 闭 | 极限点可能掉出 $K$，不动点不存在 |
> | $\rho<1$ | 例如 $K=[1,\infty)$、$\varphi(x)=x+\frac1x$ 满足 $\lvert\varphi(x)-\varphi(y)\rvert<\lvert x-y\rvert$ 却无不动点（练习 11） |
>
> 顺带得到**先验误差估计**：$\lVert x_n-z_0\rVert\le\dfrac{\rho^n}{1-\rho}\lVert x_1-x_0\rVert$ —— 这是数值迭代里判断「迭代到第几步」的标准公式。
>
> 反过来，只要把 $\rho<1$ 弱化成 $\lVert\varphi(x)-\varphi(y)\rVert<\lVert x-y\rVert$，就不再有唯一性；但**再加 $K$ 紧**，结论仍然成立（练习 11）。

### 2.4.5 练习

1. 设 $E$ 是赋范空间。证明：$E$ 是 Banach 空间 $\iff$ $E$ 中任意一列半径趋于 $0$ 的递降闭球有非空交。
2. 设 $S$ 是 $\ell^2$ 中「只有有限多个非零项」的列构成的线性子空间。证明 $S$ 不闭。
3. 设 $S$ 是非空集，$E$ 是 $\mathbb R$ 上的 Banach 空间。令 $C^b(S,E)$ 为从 $S$ 到 $E$ 的有界连续函数构成的向量空间，配范数 $\lVert f\rVert_b=\sup_{s\in S}\lVert f(s)\rVert$。证明 $C^b(S,E)$ 是 Banach 空间。
4. 设线性空间 $E$ 含一个可数集 $(x_n)_{n\in\mathbb N}$，使每个 $x\in E$ 都可写成 $x_n$ 的有限线性组合。证明：不可能给 $E$ 配上范数使 $(E,\lVert\cdot\rVert)$ 成为 Banach 空间。
5. 令

   $$
   \ell^{\infty,0}=\{(x_n)_{n\in\mathbb N}\in\ell^\infty:\ x_n\ \text{中至多只有有限多个不同的值}\}.
   $$

   问：$\ell^\infty_c$ 在 $\ell^\infty$ 中的闭包是什么？
6. 令

   $$
   S=\Bigl\{(x_n)_{n\in\mathbb N}:\ \Bigl(\sum_{k=1}^nx_k\Bigr)_{n\in\mathbb N}\ \text{收敛}\Bigr\},
   $$

   并对 $x=(x_n)_{n\in\mathbb N}\in S$ 定义 $\lVert x\rVert=\sup_{n\in\mathbb N}\bigl|\sum_{k=1}^nx_k\bigr|$。证明 $(S,\lVert\cdot\rVert)$ 是 Banach 空间。
7. 设 $E$ 是赋范空间，$Y\subseteq E$ 是闭线性子空间。若 $Y$（带 $E$ 的范数）与 $E/Y$（带商范数）都完备，证明 $E$ 完备。
8. 设 $f:C[0,1]\to\mathbb R$，$f(x)=x(1)$。
   (a) 证明若在 $C[0,1]$ 上取通常范数，则 $f$ 连续。
   (b) 证明若在 $C[0,1]$ 上取 $p$-范数 $\lVert x\rVert_p=\bigl(\int_0^1|x(t)|^p\,\mathrm dt\bigr)^{1/p}$（$1<p<\infty$），则 $f$ 不连续。
9. $BV[0,1]$ 是 Banach 空间。
10. 证明 Banach 空间 $V$ 可分 $\iff$

    $$
    V=\bigcup_{n=1}^\infty V_n,
    $$

    其中 $(V_n)_{n\in\mathbb N}$ 是 $V$ 的有限维子空间。
11. 证明：若把定理 2.4.9 的条件弱化为

    $$
    \lVert\varphi(x)-\varphi(y)\rVert<\lVert x-y\rVert,\qquad\forall x,y\in K,
    $$

    则结论不成立；但若 $K$ 紧，在此弱化条件下结论仍成立。
12. 在 $C[0,1]$ 中定义两列 $(f_n)$、$(g_n)$：

    $$
    f_n(t)=\begin{cases}1-nt,&0\le t\le\dfrac1n,\\[4pt] 0,&\dfrac1n<t\le1;\end{cases}
    \qquad
    g_n(t)=\begin{cases}n^2t,&0\le t\le\dfrac1n,\\[4pt] n(2-nt),&\dfrac1n<t\le\dfrac2n,\\[4pt] 0,&\dfrac2n<t\le1.\end{cases}
    $$

    (i) 证明 $f_n\to0$ 于 $L^1$ 范数，但 $(f_n)$ 在 $[0,1]$ 上不逐点收敛到 $0$；
    (ii) 证明 $g_n\to0$ 在 $[0,1]$ 上逐点成立，但 $\lVert g_n\rVert_{L^1}=1$，故 $(g_n)$ 不按 $L^1$ 范数收敛到 $0$。
13. 给出定理 2.4.5 的完整证明。

## 练习答案

**1.** ($\Rightarrow$) 设 $E$ 完备，$\overline B(x_n,r_n)$ 递降且 $r_n\to0$。$m>n$ 时 $x_m\in\overline B(x_n,r_n)$，故 $\lVert x_m-x_n\rVert\le r_n\to0$，$(x_n)$ 是 Cauchy 列。设 $x_n\to x$，则对每个 $n$，$x$ 是 $x_m$（$m>n$）的极限而 $\overline B(x_n,r_n)$ 闭，故 $x\in\overline B(x_n,r_n)$。

($\Leftarrow$) 设 $(x_n)$ Cauchy。取子列使

$$
\lVert x_{n_{k+1}}-x_{n_k}\rVert\le2^{-(k+1)} .
$$

则闭球 $\overline B\bigl(x_{n_k},2^{-(k-1)}\bigr)$ 递降：若 $j>k$，则

$$
\lVert x_{n_j}-x_{n_k}\rVert\le\sum_{i=k}^{j-1}2^{-(i+1)}\le2^{-k}\le2^{-(k-1)} .
$$

由假设交非空，取 $x$ 在其中。则 $\lVert x-x_{n_k}\rVert\le2^{-(k-1)}\to0$，子列收敛，故原列收敛。$\square$

**2.** 取

$$
x=\Bigl(1,\tfrac12,\tfrac13,\dots\Bigr)\in\ell^2,\qquad x^{(n)}=\Bigl(1,\tfrac12,\dots,\tfrac1n,0,0,\dots\Bigr)\in S .
$$

则 $\lVert x^{(n)}-x\rVert_2^2=\sum_{k>n}k^{-2}\to0$，而 $x\notin S$（无限多个非零项）。故 $S$ 不闭。$\square$

**3.** 设 $(f_n)$ 在 $\lVert\cdot\rVert_b$ 下 Cauchy。对每个 $s$，$\lVert f_n(s)-f_m(s)\rVert\le\lVert f_n-f_m\rVert_b\to0$，故 $(f_n(s))$ 是 $E$ 中的 Cauchy 列，$E$ 完备给出 $f(s):=\lim f_n(s)$。取 $N$ 使 $n,m>N$ 时 $\lVert f_n-f_m\rVert_b<1$，则

$$
\sup_s\lVert f(s)\rVert\le\lVert f_N\rVert_b+1<\infty,
$$

故 $f$ 有界；且令 $m\to\infty$ 得 $n>N$ 时 $\lVert f_n-f\rVert_b\le1$，即 $f_n\to f$ 于 $\lVert\cdot\rVert_b$。最后 $f$ 连续：一致极限保持连续性（若每个 $f_n$ 连续，则对 $s_\alpha\to s$，$\lVert f(s_\alpha)-f(s)\rVert\le2\lVert f-f_n\rVert_b+\lVert f_n(s_\alpha)-f_n(s)\rVert$，先取 $n$ 大使第一项小于 $\varepsilon/2$，再由 $f_n$ 连续得第二项最终小于 $\varepsilon/2$）。$\square$

> **说明** 连续性的论证只用到「一致极限保连续」，与 $S$ 的拓扑无关；若 $S$ 上取离散拓扑，所有函数都连续，本题就退化为「有界函数空间完备」。

**4.** 题面需理解为「该可数集是 $E$ 的 Hamel 基」（否则 $E=\mathbb K$、$\{x_n\}\equiv1$ 是反例，而 $\mathbb K$ 当然是 Banach 空间）。此时 $E$ 无限维，且

$$
E=\bigcup_{n=1}^\infty F_n,\qquad F_n=\operatorname{span}\{x_1,\dots,x_n\}.
$$

每个 $F_n$ 是有限维子空间，因而**闭**（命题 1.5.19）；又 $F_n\ne E$（Hamel 基无限），故 $F_n$ 是闭的真子空间，内部为空（若含球则含 $0$ 的邻域，由 2.1.9 练习 11 得 $F_n=E$）。于是 $E$ 是可数个无处稠密集之并。若 $E$ 完备，Baire 纲定理要求其中某个 $F_n$ 有内点，矛盾。故 $E$ 不能是 Banach 空间。$\square$

**5.** 题面有两处记号不一致（定义 $\ell^{\infty,0}$、提问 $\ell^\infty_c$），两种读法答案不同：

**(一) 若 $\ell^\infty_c$ 指 $c_{00}$（有限支撑）**，则 $\overline{c_{00}}=c_0$。$\subseteq$：$c_{00}\subseteq c_0$ 且 $c_0$ 闭。$\supseteq$：$x\in c_0$ 时截断 $x^{(n)}=(\xi_1,\dots,\xi_n,0,\dots)\in c_{00}$ 且 $\lVert x-x^{(n)}\rVert_\infty=\sup_{k>n}\lvert\xi_k\rvert\to0$。

**(二) 若按字面「只有有限多个不同的取值」**，则该集合 $A$ **在 $\ell^\infty$ 中稠密**，故闭包是整个 $\ell^\infty$。证明：任取 $x\in\ell^\infty$ 与 $\varepsilon>0$，因 $\lVert x\rVert_\infty\le M$，把每个坐标四舍五入到网格 $\varepsilon\mathbb Z\cap[-M-\varepsilon,M+\varepsilon]$，得 $r(x)$ 只取**有限多个**值，且 $\lVert x-r(x)\rVert_\infty\le\varepsilon$。$\square$

> **要点** 两种读法的反差正说明记号必须写清：$c_{00}$ 的闭包是「小」的 $c_0$，而「有限多值」的闭包是「大」的 $\ell^\infty$。教材把两者混成了一个符号。

**6.** 令 $\Phi(x)=\bigl(\sum_{k\le n}x_k\bigr)_{n\in\mathbb N}$。**值域**：$\Phi(x)\in c$ 当且仅当部分和收敛，正是 $S$ 的定义，故 $\Phi:S\to c$ 是线性双射（$x_n=s_n-s_{n-1}$ 给出逆）。**等距**：$\lVert\Phi(x)\rVert_\infty=\sup_n\lvert s_n\rvert=\lVert x\rVert$。$c$ 是 Banach 空间，故 $S$ 也是。$\square$

**7.** 设 $(x_n)$ 是 $E$ 中 Cauchy 列。则 $([x_n])$ 是 $E/Y$ 中 Cauchy 列，设 $[x_n]\to[x]$。对每个 $n$ 取 $u_n\in Y$ 使

$$
\lVert x_n-x-u_n\rVert\le\lVert[x_n]-[x]\rVert+\frac1n .
$$

则 $(u_n)$ 在 $Y$ 中 Cauchy：

$$
\lVert u_n-u_m\rVert\le\lVert u_n-(x_n-x)\rVert+\lVert x_n-x_m\rVert+\lVert(x_m-x)-u_m\rVert\longrightarrow0 .
$$

$Y$ 完备，设 $u_n\to u\in Y$。于是

$$
x_n=(x_n-x-u_n)+u_n+x\longrightarrow 0+u+x\in E .
$$

故 $E$ 完备。$\square$

**8.** **(a)** $\lvert f(x)\rvert=\lvert x(1)\rvert\le\lVert x\rVert_\infty$，故 $\lVert f\rVert\le1$（取 $x\equiv1$ 得 $\lVert f\rVert=1$）。

**(b)** 取 $x_n(t)=t^n$，则 $\lVert x_n\rVert_p^p=\int_0^1t^{np}dt=\dfrac{1}{np+1}\to0$，即 $x_n\to0$；但 $f(x_n)=1\not\to0=f(0)$。故 $f$ 不连续。$\square$

> **要点** 同一个线性泛函在两个范数下连续性相反，说明 $C[0,1]$ 上 $\lVert\cdot\rVert_p$ 与 $\lVert\cdot\rVert_\infty$ **不等价**（且 $\lVert\cdot\rVert_p$ 下 $C[0,1]$ 甚至不完备，其完备化是 $L^p[0,1]$）。

**9.** 记

$$
\lVert x\rVert_{BV}=\lvert x(0)\rvert+V(x),\qquad V(x)=\sup\Bigl\{\sum_{i=1}^n\lvert x(t_i)-x(t_{i-1})\rvert:0=t_0<\dots<t_n=1\Bigr\}.
$$

**是范数。** 正定性与齐次性显然；三角不等式由 $V(x+y)\le V(x)+V(y)$（逐项三角不等式后取上确界）。

**完备。** 设 $(x_k)$ 在 $\lVert\cdot\rVert_{BV}$ 下 Cauchy。首先

$$
\lvert x_k(t)-x_j(t)\rvert\le\lvert x_k(0)-x_j(0)\rvert+\lvert(x_k-x_j)(t)-(x_k-x_j)(0)\rvert\le\lVert x_k-x_j\rVert_{BV},
$$

故每个 $t$ 处 $(x_k(t))$ 是 Cauchy 列，设 $x_k(t)\to x(t)$。其次，对任意分点组 $\{t_i\}$，

$$
\sum_i\lvert x(t_i)-x(t_{i-1})\rvert=\lim_k\sum_i\lvert x_k(t_i)-x_k(t_{i-1})\rvert\le\sup_kV(x_k)<\infty,
$$

故 $V(x)\le\sup_kV(x_k)<\infty$，$x\in BV$。最后，对 $\varepsilon>0$ 取 $K$ 使 $k,j>K$ 时 $\lVert x_k-x_j\rVert_{BV}<\varepsilon$，对任意分点组令 $j\to\infty$ 得 $V(x_k-x)\le\varepsilon$，再结合 $t=0$ 处的收敛得 $\lVert x_k-x\rVert_{BV}\le\varepsilon$。故 $x_k\to x$。$\square$

**10.** **题面按字面是错的。** 若 $V=\bigcup_nV_n$（真并），则 $V$ 是可数个闭真子空间之并；每个 $V_n$ 无处稠密，与 Baire 纲定理矛盾（$V$ 完备且无限维时）。正确的说法是把等号换成**闭包**：

$$
V\ \text{可分}\iff V=\overline{\bigcup_{n=1}^\infty V_n}\quad(V_n\ \text{有限维}).
$$

**($\Leftarrow$)** 每个有限维 $V_n$ 可分（取有理系数组合即可数稠密），记 $D_n$ 为其可数稠密集。则 $D=\bigcup_nD_n$ 可数，且对 $x\in V$、$\varepsilon>0$，由 $\bigcup V_n$ 稠密取 $y\in V_n$ 使 $\lVert x-y\rVert<\varepsilon/2$，再取 $d\in D_n$ 使 $\lVert y-d\rVert<\varepsilon/2$，故 $D$ 稠密。

**($\Rightarrow$)** 设 $D=\{d_n\}$ 可数稠密，令 $V_n=\operatorname{span}\{d_1,\dots,d_n\}$。则 $\bigcup_nV_n=\operatorname{span}D$，而 $\overline{\operatorname{span}}D=V$（$D$ 稠密），故 $V=\overline{\bigcup_nV_n}$。$\square$

**11.** **弱化后一般不再有不动点。** 取 $K=[1,\infty)\subseteq\mathbb R$（闭、无界）、$\varphi(x)=x+\dfrac1x$。则 $\varphi(K)\subseteq K$，且由中值定理

$$
\lvert\varphi(x)-\varphi(y)\rvert=\Bigl\lvert1-\frac1{\xi^2}\Bigr\rvert\lvert x-y\rvert<\lvert x-y\rvert\qquad(x\ne y),
$$

但 $\varphi(x)>x$，无不动点。

**$K$ 紧时仍有不动点。** 令 $g(x)=\lVert\varphi(x)-x\rVert$，它连续，$K$ 紧故在 $x_0$ 处取到最小值 $m$。若 $m>0$，则由假设

$$
g\bigl(\varphi(x_0)\bigr)=\lVert\varphi(\varphi(x_0))-\varphi(x_0)\rVert<\lVert\varphi(x_0)-x_0\rVert=g(x_0)=m,
$$

与 $m$ 的最小性矛盾（注意 $\varphi(x_0)\in K$）。故 $m=0$，$x_0$ 是不动点。**唯一性也保留**：若 $u\ne v$ 都是不动点，则 $\lVert u-v\rVert=\lVert\varphi(u)-\varphi(v)\rVert<\lVert u-v\rVert$，矛盾。$\square$

**12.** **(i)** $f_n$ 的图形是底 $\frac1n$、高 $1$ 的三角形，故

$$
\lVert f_n\rVert_{L^1}=\frac12\cdot\frac1n\cdot1=\frac{1}{2n}\longrightarrow0 .
$$

但 $f_n(0)=1$ 对一切 $n$ 成立，故 $(f_n(0))$ 不收敛到 $0$。

**(ii)** 固定 $t\in(0,1]$，当 $n>2/t$ 时 $2/n<t$，此时 $g_n(t)=0$；而 $g_n(0)=0$。故 $g_n(t)\to0$ 对每个 $t\in[0,1]$ 成立。但 $g_n$ 的图形是底 $\frac2n$、高 $n$ 的三角形，

$$
\lVert g_n\rVert_{L^1}=\frac12\cdot\frac2n\cdot n=1\not\longrightarrow0 .
$$

$\square$

> **要点** 两条方向都断：$L^1$ 收敛不推逐点收敛（(i)）；逐点收敛不推 $L^1$ 收敛（(ii)，缺控制函数）。这正是控制收敛定理（Dominated Convergence）里「被可积函数控制」不可省的原因。

**13.** **构造。** 设 $\mathcal C$ 是 $E$ 中全体 Cauchy 列的集合（逐项运算使其成为向量空间），$\mathcal N\subseteq\mathcal C$ 是零列（$x_n\to0$）构成子空间。对 $(x_n)\in\mathcal C$，$(\lVert x_n\rVert)$ 是 $\mathbb R$ 中的 Cauchy 列，故

$$
p\bigl((x_n)\bigr):=\lim_{n\to\infty}\lVert x_n\rVert
$$

存在；$p$ 是 $\mathcal C$ 上的半范数，且 $p^{-1}(0)=\mathcal N$。令

$$
\widetilde E=\mathcal C/\mathcal N,\qquad \bigl\lVert[(x_n)]\bigr\rVert:=p\bigl((x_n)\bigr).
$$

这是范数（半范数商掉核）。

**$T$ 是等距。** $T(x)=[(x,x,x,\dots)]$，则 $\lVert Tx\rVert=\lVert x\rVert$。

**像稠密。** 设 $\xi=[(x_n)]\in\widetilde E$。由 Cauchy 性，取子列使 $\lVert x_m-x_{n_k}\rVert<2^{-k}$（$m>n_k$），则

$$
\lVert\xi-Tx_{n_k}\rVert=\lim_{m}\lVert x_m-x_{n_k}\rVert\le2^{-k}\longrightarrow0 .
$$

**完备。** 设 $(\xi_k)$ 在 $\widetilde E$ 中 Cauchy。取子列使 $\lVert\xi_{k+1}-\xi_k\rVert<2^{-(k+1)}$，再由商范数的定义逐次挑选代表元使

$$
\lVert x^{(k+1)}-x^{(k)}\rVert<2^{-k}
$$

（先取 $x^{(1)}$，已知 $x^{(k)}$ 时取 $x^{(k+1)}\in\xi_{k+1}$ 使 $\lVert x^{(k+1)}-x^{(k)}\rVert<\lVert\xi_{k+1}-\xi_k\rVert+2^{-(k+1)}<2^{-k}$）。则 $(x^{(k)})$ 在 $E$ 中 Cauchy：$k<l$ 时

$$
\lVert x^{(l)}-x^{(k)}\rVert\le\sum_{i=k}^{l-1}2^{-i}\le2^{-(k-1)} .
$$

记 $\xi=[(x^{(k)})]$，则 $\lVert\xi_k-\xi\rVert\le2^{-(k-1)}\to0$，故 $\widetilde E$ 完备。于是 $\widetilde E$ 是 Banach 空间、$T$ 是等距、$T(E)$ 稠密。

**唯一性。** 设 $\widetilde E_1,\widetilde E_2$ 都是完备化，$T_i:E\to\widetilde E_i$。在 $T_1(E)$ 上令 $S(T_1x)=T_2x$，则 $S$ 是等距，等距当然一致连续，而 $\widetilde E_2$ 完备，故 $S$ 唯一延拓为 $\overline{T_1(E)}=\widetilde E_1$ 上的等距；其像含稠密集 $T_2(E)$ 且闭（等距像），故为满射。于是 $\widetilde E_1\cong\widetilde E_2$。$\square$

> **要点** 这个构造把「不完备」理解成「有 Cauchy 列没极限」，于是干脆**把 Cauchy 列本身当作新点**，再把互相靠拢的列（差为零列）粘成一点。有理数 $\mathbb Q$ 到 $\mathbb R$ 的 Cantor 构造是它的特例。

---

## 订正

1. **定理 2.4.9 证明末句的映射写错。** 原文在导出 $z_0=\varphi(z_0)$ 时写「$\varphi(x_n)\to\varphi(x_0)$」；$x_n\to z_0$ 才是已证事实，故应为 $\varphi(x_n)\to\varphi(z_0)$。（若真取 $x_0$，只能得到 $z_0=\varphi(x_0)$，与 $x_0$ 的任意性矛盾。）
2. **定理 2.4.9 的迭代起点与递推式不一致。** 原文写「set $x_{n+1}=\varphi(x_n)$ for $n=1,2,\dots$」，但其后立即用到 $x_0$ 与 $\lVert x_1-x_0\rVert$，故递推应从 $n=0$ 起。本文已改。
3. **定理 2.4.4 证明中那个级数下标从 $1$ 起是笔误。** 原文写 $x_{n_1}+\sum_{k=1}^\infty(x_{n_k}-x_{n_{k-1}})$，但 $x_{n_0}$ 未定义；应为 $\sum_{k=2}^\infty$。本文已改。
4. **定理 2.4.5 的「唯一」未给证明**（原文指向参考书）。本文在 2.4.2 补充里给出标准论证：在稠密子集上定义等距，再用完备性延拓。
5. **练习 5 的记号自相矛盾**：定义的是 $\ell^{\infty,0}$，提问却写 $\ell^\infty_c$；且按字面「只有有限多个不同的取值」定义的集合不是线性子空间。本文在答案里把两种读法都算了：按 $c_{00}$ 理解闭包是 $c_0$，按字面理解闭包是整个 $\ell^\infty$。
6. **练习 10 的等号应为闭包。** 「$V=\bigcup_{n=1}^\infty V_n$」在 $V$ 无限维时**不可能**成立：有限维子空间闭且无处稠密，可数并仍无处稠密，与 $V$ 完备（Baire）矛盾。正确说法是 $V=\overline{\bigcup_nV_n}$，本文已按此证明。
7. **练习 4 缺一个前提。** 若 $E$ 有限维，题面的可数集可以取成 $\{x_n\}\equiv x_0\ne0$，此时 $E$ 当然是 Banach 空间 —— 题面应理解为「该可数集是 $E$ 的 Hamel 基」。

---

## 下一步

第 2 章至此结束。第 3 章起进入**另一半支柱：Hilbert 空间** —— 给赋范空间再装一个内积，于是「长度」之外还有「夹角」，最佳逼近、正交分解、正交基全都随之而来。

---

系列：[[01-预备·Zorn引理与向量空间]] · [[02-度量空间与拓扑空间]] · [[03-拓扑向量空间]] · [[04-赋范空间]] · [[05-Banach空间]] · [[06-可分Banach空间]] · [[08-Hilbert空间·基本性质]] · [[09-Hilbert空间·最佳逼近]] · [[10-Hilbert空间·正交分解]] · [[11-Hilbert空间·正交基]] · [[12-连续线性算子·有界性与算子范数]] · [[13-连续线性算子·连续线性泛函与弱拓扑]] · [[14-连续线性算子·对偶空间的例子]] · [[15-连续线性算子·伴随算子]] · [[16-连续线性算子·投影算子]]
