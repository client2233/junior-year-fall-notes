# 08 Hilbert 空间 · 基本性质

> 系列：[[01-预备·Zorn引理与向量空间]] · [[02-度量空间与拓扑空间]] · [[03-拓扑向量空间]] · [[04-赋范空间]] · [[05-Banach空间]] · [[06-可分Banach空间]] · [[07-Banach空间·完备化与紧性]] · [[09-Hilbert空间·最佳逼近]] · [[10-Hilbert空间·正交分解]] · [[11-Hilbert空间·正交基]] · [[12-连续线性算子·有界性与算子范数]] · [[13-连续线性算子·连续线性泛函与弱拓扑]] · [[14-连续线性算子·对偶空间的例子]] · [[15-连续线性算子·伴随算子]] · [[16-连续线性算子·投影算子]]
>
> 第 2 章把「长度」抽象成了**范数**。第 3 章再往上抽象一层：$\mathbb R^n$ 之所以好用，不只是因为向量有长度，还因为有**夹角**。夹角来自**内积**，而内积同时给出长度（$\lVert x\rVert^2=\langle x,x\rangle$）。所以内积比范数**更基本** —— 这就是 Hilbert 空间比 Banach 空间更接近 $\mathbb R^n$ 的原因。
>
> 两条对照，先记住：
>
> - 「坐标」这一半：Schauder 基把 $\mathbb R^n$ 的坐标推广到（部分）可分 Banach 空间，但并非每个可分 Banach 空间都有 Schauder 基（[[06-可分Banach空间]]）。而 **Hilbert 空间一定有基**（正交基，见 3.4）。
> - 「长度」这一半：范数由内积导出有**判据**（平行四边形律，定理 3.1.10），并非任意范数都行。
>
> **记号约定**：本书内积取「对第一变量线性、对第二变量共轭线性」的约定，即 $\langle\lambda x,y\rangle=\lambda\langle x,y\rangle$、$\langle x,\lambda y\rangle=\overline\lambda\langle x,y\rangle$。相应地 $\mathbb C^n$、$L^2$、$\ell^2$ 的内积共轭落在**第二个**变量上。原 PDF 的文字层不携带共轭上划线，本文一律按此约定补齐。

## 3.1 Hilbert 空间的基本性质

### 3.1.1 内积的定义

**定义 3.1.1（内积）** 设 $E$ 是 $\mathbb K$ 上的向量空间。函数

$$
\langle\cdot,\cdot\rangle:E\times E\to\mathbb K
$$

称 $E$ 上的**内积**，若

1. $\langle x,x\rangle\ge0$ 对一切 $x\in E$；且 $\langle x,x\rangle=0$ 当且仅当 $x=0$；
2. $\langle x,y\rangle=\overline{\langle y,x\rangle}$ 对一切 $x,y\in E$；
3. $\langle\lambda x,y\rangle=\lambda\langle x,y\rangle$ 对一切 $x,y\in E$；
4. $\langle x+y,z\rangle=\langle x,z\rangle+\langle y,z\rangle$ 对一切 $x,y,z\in E$。

当底域 $\mathbb K=\mathbb R$ 或 $\mathbb C$ 时，分别称 $(E,\langle\cdot,\cdot\rangle)$ 为**实内积空间**或**复内积空间**。$\langle x,y\rangle$ 称 $x$ 与 $y$ 的**内积**。

> **说明（为什么共轭落在第二个变量）** 条件 (2) 与 (3) 一起决定了内积只能「一头线性、一头共轭线性」：在复空间里，$\langle\lambda x,\lambda x\rangle=\lambda\overline{\lambda}\langle x,x\rangle=|\lambda|^2\lVert x\rVert^2$ 才能与 (1) 相容。至于共轭放在哪一头，是**约定**：本书放在第二个变量上，因此 (3) 无共轭、(3') 有共轭。另一套流行的约定（物理文献）恰好相反，引用结论时须先对表。

**注 3.1.2** 内积空间 $(E,\langle\cdot,\cdot\rangle)$ 还自然满足

3'. $\langle x,\lambda y\rangle=\overline\lambda\langle x,y\rangle$ 对一切 $x,y\in E$ 与 $\lambda\in\mathbb K$；
4'. $\langle z,x+y\rangle=\langle z,x\rangle+\langle z,y\rangle$ 对一切 $x,y,z\in E$。

**证明** 由 (2) 与 (3)：

$$
\langle x,\lambda y\rangle=\overline{\langle\lambda y,x\rangle}=\overline{\lambda\langle y,x\rangle}=\overline\lambda\,\overline{\langle y,x\rangle}=\overline\lambda\langle x,y\rangle,
$$

同法可得 (4')。$\square$

> **要点** (3) 与 (3') 的差别只有一条**共轭横线**，这是复内积空间里最容易出错的地方。检验自己是否记牢：$\langle ix,y\rangle=i\langle x,y\rangle$，但 $\langle x,iy\rangle=\overline i\langle x,y\rangle=-i\langle x,y\rangle$。

### 3.1.2 Cauchy–Schwarz 不等式

**定理 3.1.3（Cauchy–Schwarz 不等式）** 设 $(E,\langle\cdot,\cdot\rangle)$ 是内积空间，则

$$
|\langle x,y\rangle|\le\sqrt{\langle x,x\rangle}\,\sqrt{\langle y,y\rangle},\qquad\forall x,y\in E.
$$

等号成立当且仅当存在 $\lambda\in\mathbb K$ 使 $x=\lambda y$，或 $y=0$。

**证明** 若 $y=0$ 等号显然。设 $y\ne0$。对任意 $\mu\in\mathbb K$，

$$
0\le\langle x-\mu y,x-\mu y\rangle=\langle x,x\rangle-\mu\langle y,x\rangle-\overline\mu\langle x,y\rangle+|\mu|^2\langle y,y\rangle.
$$

取 $\mu=\dfrac{\langle x,y\rangle}{\langle y,y\rangle}$，代入得

$$
0\le\langle x,x\rangle-\frac{|\langle x,y\rangle|^2}{\langle y,y\rangle}-\frac{|\langle x,y\rangle|^2}{\langle y,y\rangle}+\frac{|\langle x,y\rangle|^2}{\langle y,y\rangle}
=\langle x,x\rangle-\frac{|\langle x,y\rangle|^2}{\langle y,y\rangle},
$$

移项即得所需不等式。

当 $x=\lambda y$ 时等号显然。反之设等号成立，则对上述 $\mu$ 有

$$
\langle x-\mu y,x-\mu y\rangle=0,
$$

由定义 (1) 得 $x=\mu y$。$\square$

> **补充（等号条件的两种情形）** 定理陈述写成「$x=\lambda y$ 或 $y=0$」，其实两种情形可以统一：等号成立 $\iff$ $x$ 与 $y$ **线性相关**（$y=0$ 时任何 $x$ 都与 $y$ 相关）。分开写只是为了回避「$y=0$ 时 $\mu$ 没有定义」这个技术细节。
>
> **补充（Cauchy–Schwarz 直接给出三角不等式）** 由 $\lVert x\rVert^2=\langle x,x\rangle$ 与 Cauchy–Schwarz，

$$
\lVert x+y\rVert^2=\lVert x\rVert^2+2\operatorname{Re}\langle x,y\rangle+\lVert y\rVert^2\le\lVert x\rVert^2+2\lVert x\rVert\lVert y\rVert+\lVert y\rVert^2=(\lVert x\rVert+\lVert y\rVert)^2,
$$

> 故 $\lVert x+y\rVert\le\lVert x\rVert+\lVert y\rVert$。这是 Cauchy–Schwarz 在所有分析分支中最重要的一次应用，也说明「由内积导出的 $\lVert\cdot\rVert$ 是范数」中唯一不平凡的一条就是三角不等式。

### 3.1.3 由内积导出的范数

**定义 3.1.4（诱导范数与 Hilbert 空间）** 设 $\langle\cdot,\cdot\rangle$ 是向量空间 $E$ 上的内积。定义 $E$ 上的范数

$$
\lVert x\rVert=\sqrt{\langle x,x\rangle},\qquad\forall x\in E.
$$

$E$ 的范数拓扑也称 $E$ 的**内积拓扑**。若 $E$ 在上述范数下完备，则称 $E$ 为 **Hilbert 空间**。

> **补充（三角不等式要补一步）** 教材说「易验证 $\lVert\cdot\rVert$ 满足定义 2.1.1 的范数条件」。其中 (1) 正定性直接来自内积定义 (1)；(2) 齐性由 $\lVert\lambda x\rVert^2=|\lambda|^2\lVert x\rVert^2$ 得到；(3) 三角不等式**并不显然**，它就是 3.1.2 补充里那一步 —— 必须先有 Cauchy–Schwarz。本文已把它补齐。

### 3.1.4 标准例子

**例 3.1.5（$\mathbb R^n$ 与 $\mathbb C^n$）** $\mathbb R^n$ 上的内积定义为

$$
\langle x,y\rangle=\sum_{i=1}^nx_iy_i,\qquad\forall x=(x_1,\dots,x_n)^t,\ y=(y_1,\dots,y_n)^t\in\mathbb R^n;
$$

$\mathbb C^n$ 上的内积定义为

$$
\langle x,y\rangle=\sum_{i=1}^nx_i\overline{y_i},\qquad\forall x=(x_1,\dots,x_n)^t,\ y=(y_1,\dots,y_n)^t\in\mathbb C^n.
$$

**例 3.1.6（$L^2[a,b]$）**

$$
\langle x,y\rangle=\int_a^bx(t)\overline{y(t)}\,\mathrm dt,\qquad\forall x,y\in L^2[a,b].
$$

**例 3.1.7（$\ell^2$）**

$$
\langle x,y\rangle=\sum_{n=1}^\infty x_n\overline{y_n},\qquad\forall x=(x_n)_{n\in\mathbb N},\ y=(y_n)_{n\in\mathbb N}\in\ell^2.
$$

**例 3.1.8（矩阵空间 $M_{m,n}$）**

$$
\langle A,B\rangle=\operatorname{tr}(B^*A),\qquad\forall A,B\in M_{m,n},
$$

其中 $B^*=(\overline{b_{ji}})$ 是 $B=(b_{ij})\in M_{m,n}$ 的**共轭转置**（伴随矩阵），$\operatorname{tr}(\cdot)$ 是 $M_n:=M_{n,n}$ 上的**迹泛函**

$$
\operatorname{tr}(C)=\sum_{i=1}^nc_{ii},\qquad\forall C=(c_{ij})_{1\le i,j\le n}\in M_n.
$$

它诱导的正是 [[04-赋范空间]] 例 2.1.4 中的 **Frobenius 范数** $\lVert(a_{ij})\rVert=\bigl(\sum_{i,j}|a_{ij}|^2\bigr)^{1/2}$。

> **补充（$\operatorname{tr}(B^*A)$ 为什么等于逐项内积）** 直接算：

$$
\operatorname{tr}(B^*A)=\sum_{i=1}^n(B^*A)_{ii}=\sum_{i=1}^n\sum_{k=1}^m(B^*)_{ik}A_{ki}=\sum_{i,k}\overline{b_{ki}}a_{ki},
$$

> 即把 $A$、$B$ 拉成 $mn$ 维向量后的标准内积。所以 $M_{m,n}$ 与 $\mathbb K^{mn}$ 是**等距同构**的内积空间，Frobenius 范数就是 $\mathbb K^{mn}$ 上的 $2$-范数。

### 3.1.5 内积的连续性

**定理 3.1.9（内积的连续性）** 设 $x_n,y_n,x,y$ 是内积空间 $E$ 中的元素（$n=1,2,\dots$），且在内积拓扑下 $x_n\to x$、$y_n\to y$。则

$$
\langle x_n,y_n\rangle\to\langle x,y\rangle.
$$

**证明** 因 $y_n\to y$，$(y_n)$ 有界，即存在 $C>0$ 使 $\lVert y_n\rVert\le C$（$n\in\mathbb N$）。由 Cauchy–Schwarz 与三角不等式，

$$
|\langle x,y\rangle-\langle x_n,y_n\rangle|\le|\langle x,y\rangle-\langle x,y_n\rangle|+|\langle x,y_n\rangle-\langle x_n,y_n\rangle|\le\lVert x\rVert\lVert y-y_n\rVert+\lVert x-x_n\rVert\lVert y_n\rVert\le\lVert x\rVert\lVert y-y_n\rVert+C\lVert x-x_n\rVert,
$$

右端随 $n\to\infty$ 趋于 $0$，故 $|\langle x,y\rangle-\langle x_n,y_n\rangle|\to0$。$\square$

> **要点** 这条定理的用处是「**极限可以搬进内积**」：只要两个变量都收敛，内积就收敛，无需任何附加条件。后面正交分解、正交基、投影算子的构造全靠它。
>
> 注意它**对两个变量同时成立**（不是逐个变量），证明里用 $\lVert y_n\rVert$ 有界来把双变量误差拆成两份单变量误差。

### 3.1.6 平行四边形律与极化恒等式

**定理 3.1.10（平行四边形律）** 向量空间 $E$ 上的范数 $\lVert\cdot\rVert$ 由 $E$ 上的某个内积导出，当且仅当

$$
\lVert x+y\rVert^2+\lVert x-y\rVert^2=2(\lVert x\rVert^2+\lVert y\rVert^2),\qquad\forall x,y\in E. \tag{3.1}
$$

**证明**

**（$\Rightarrow$）** 若 $\lVert\cdot\rVert$ 由内积导出，直接展开：

$$
\lVert x\pm y\rVert^2=\lVert x\rVert^2\pm\langle x,y\rangle\pm\langle y,x\rangle+\lVert y\rVert^2,
$$

两式相加即得 (3.1)。

**（$\Leftarrow$）** 设 (3.1) 对一切 $x,y\in E$ 成立，要**造出**那个内积。假想 $E$ 已经是内积空间，先做观察。

**实情形。** 由展开式

$$
\lVert x+y\rVert^2=\lVert x\rVert^2+\langle x,y\rangle+\langle y,x\rangle+\lVert y\rVert^2,\qquad
\lVert x-y\rVert^2=\lVert x\rVert^2-\langle x,y\rangle-\langle y,x\rangle+\lVert y\rVert^2,
$$

相减得 $\lVert x+y\rVert^2-\lVert x-y\rVert^2=2\langle x,y\rangle+2\langle y,x\rangle=4\langle x,y\rangle$（实情形 $\langle x,y\rangle=\langle y,x\rangle$），即

$$
\langle x,y\rangle=\frac14\bigl(\lVert x+y\rVert^2-\lVert x-y\rVert^2\bigr). \tag{3.2}
$$

**复情形。** 对 $k=0,1,2,3$，

$$
\lVert x+i^ky\rVert^2=\langle x+i^ky,x+i^ky\rangle=\lVert x\rVert^2+i^k\langle y,x\rangle+i^{-k}\langle x,y\rangle+\lVert y\rVert^2.
$$

于是

$$
\begin{aligned}
\sum_{k=0}^3i^k\lVert x+i^ky\rVert^2&=(\underbrace{1+i-1-i}_{=0})\lVert x\rVert^2+(\underbrace{1-1+1-1}_{=0})\langle y,x\rangle\\
&\qquad+(\underbrace{1+1+1+1}_{=4})\langle x,y\rangle+(\underbrace{1+i-1-i}_{=0})\lVert y\rVert^2=4\langle x,y\rangle,
\end{aligned}
$$

即

$$
\langle x,y\rangle=\frac14\sum_{k=0}^3i^k\lVert x+i^ky\rVert^2. \tag{3.3}
$$

(3.2) 与 (3.3) 称为**极化恒等式**（polarization identity）。下面就用它来**定义** $E$ 上的内积，并验证它确实是内积 —— 这一步要用到平行四边形律 (3.1)。只做复情形，实情形更简单。

**第 1 步：先分离实部与虚部。** 由 (3.3) 的实部（注意 $k=0,3$ 两项的系数为实、$k=1,2$ 两项相互抵消的虚部）得到

$$
\operatorname{Re}\langle x,y\rangle=\frac14\bigl(\lVert x+y\rVert^2-\lVert x-y\rVert^2\bigr), \tag{3.4}
$$

又 $\lVert ix\pm y\rVert=\lVert x\mp iy\rVert$（因为 $ix+y=i(x-iy)$、$ix-y=i(x+iy)$，而 $\lVert i(\cdot)\rVert=\lVert\cdot\rVert$），故

$$
\langle x,y\rangle=\operatorname{Re}\langle x,y\rangle-i\operatorname{Re}\langle ix,y\rangle. \tag{3.5}
$$

**第 2 步：可加性。** 由 (3.4) 与 (3.1)：

$$
\begin{aligned}
\operatorname{Re}\langle x+z,y\rangle+\operatorname{Re}\langle x-z,y\rangle
&=\frac14\bigl(\lVert x+z+y\rVert^2-\lVert x+z-y\rVert^2+\lVert x-z+y\rVert^2-\lVert x-z-y\rVert^2\bigr)\\
&=\frac12\bigl(\lVert x+y\rVert^2+\lVert z\rVert^2-\lVert x-y\rVert^2-\lVert z\rVert^2\bigr)=2\operatorname{Re}\langle x,y\rangle.
\end{aligned}\tag{3.6}
$$

（第二个等号：对 $\lVert(x+y)\pm z\rVert^2$ 与 $\lVert(x-y)\pm z\rVert^2$ 各用一次平行四边形律再合并。）

在 (3.6) 中取 $z=x$，得

$$
\operatorname{Re}\langle 2x,y\rangle=2\operatorname{Re}\langle x,y\rangle,
$$

于是 (3.6) 可改写为

$$
\operatorname{Re}\langle x+z,y\rangle+\operatorname{Re}\langle x-z,y\rangle=\operatorname{Re}\langle 2x,y\rangle. \tag{3.7}
$$

在 (3.7) 中把 $x$ 换成 $\dfrac12(x+z)$、$z$ 换成 $\dfrac12(x-z)$，得

$$
\operatorname{Re}\langle x,y\rangle+\operatorname{Re}\langle z,y\rangle=\operatorname{Re}\langle x+z,y\rangle. \tag{3.8}
$$

再把 (3.8) 中的 $x,z$ 换成 $ix,iz$，并用 (3.5)，得

$$
\langle x,y\rangle+\langle z,y\rangle=\langle x+z,y\rangle.
$$

**第 3 步：齐性。** 上式立即给出「正整数倍」的齐性，逐次可推到 $\alpha\in\mathbb Q$。再由范数函数连续（$q_n\to\alpha\Rightarrow\lVert(q_n-\alpha)x\rVert\to0$，配合 (3.4) 中 $\lVert(q_n-\alpha)x\pm y\rVert^2\to\lVert\pm y\rVert^2$），$x\mapsto\langle x,y\rangle$ 关于纯量因子连续，故 $\langle\alpha x,y\rangle=\alpha\langle x,y\rangle$ 对一切 $\alpha\in\mathbb R$ 成立。又由构造 $\langle ix,y\rangle=i\langle x,y\rangle$ 是直接的，故

$$
\langle\alpha x,y\rangle=\alpha\langle x,y\rangle,\qquad\forall\alpha\in\mathbb C .
$$

**第 4 步：其余公理。** 由 (3.4) 显然 $\langle x,x\rangle=\lVert x\rVert^2\ge0$，且 $\langle x,x\rangle=0\iff x=0$。而 $\langle x,y\rangle=\overline{\langle y,x\rangle}$：实部由 (3.4) 对称，虚部由 $\operatorname{Re}\langle iy,x\rangle$ 与 $-\operatorname{Re}\langle ix,y\rangle$ 相等得到（两者都等于 $\operatorname{Im}\langle x,y\rangle$）。$\square$

> **要点**
>
> 1. **定理把「内积」换成了「一条可以逐点验证的恒等式」。** 判断一个范数是否来自内积，只需在 $x,y$ 上试 (3.1)，不必去找内积。这条定理称 **Jordan–von Neumann 定理**（与极化恒等式合起来），是 Hilbert 空间理论的入口。
> 2. **正交性也是内积的产物**：$\operatorname{Re}\langle x,y\rangle=0$ 加上虚部就完全决定 $\langle\cdot,\cdot\rangle$，所以 (3.2)–(3.5) 说明「一个范数落在 $2$ 维子空间上的形状」就足以确定整个内积。
> 3. **注意反例的用处**：由练习 (1)(2)，$\lVert\cdot\rVert_\infty$ 与 $\lVert\cdot\rVert_p$（$p\ne2$）都**不满足** (3.1)，故 $C[0,1]$、$L^p$（$p\ne2$）都不是内积空间 —— 这是 $L^2$ 与 $\ell^2$ 在所有 $L^p$、$\ell^p$ 中「独一份」的根源。

### 3.1.7 等距与内积保持

**命题 3.1.11** 设 $T:H\to K$ 是内积空间之间的**线性**映射。则 $T$ 是等距当且仅当 $T$ 保持内积，即

$$
\lVert Tx\rVert=\lVert x\rVert\ \ (\forall x\in H)\iff\langle Tx,Ty\rangle=\langle x,y\rangle\ \ (\forall x,y\in H).
$$

注意两个条件都**不**保证 $T$ 是满射。

**证明** 若 $T$ 保持内积，则 $\lVert Tx\rVert^2=\langle Tx,Tx\rangle=\langle x,x\rangle=\lVert x\rVert^2$，$T$ 是等距。

反之设 $T$ 保持范数。由极化恒等式 (3.2)(3.3)（$T$ 线性，故 $T(x+\lambda y)=Tx+\lambda\,Ty$）：

$$
\langle Tx,Ty\rangle=\frac14\sum_{k=0}^3i^k\lVert Tx+i^kTy\rVert^2=\frac14\sum_{k=0}^3i^k\lVert x+i^ky\rVert^2=\langle x,y\rangle .
$$

（实情形用 (3.2)。）$\square$

> **要点** 「保范 $\Rightarrow$ 保内积」是**线性**映射独有的福利（极化恒等式要用 $T(x+y)=Tx+Ty$）。非线性等距（如平移）保范但不保内积。这也说明**引入内积并没有增加公理负担**：在线性映射这一层，长度信息已经内含了夹角信息。

### 3.1.8 内积空间的直和

**命题 3.1.12** 设 $H$、$K$ 是内积空间，内积分别为 $\langle\cdot,\cdot\rangle_H$、$\langle\cdot,\cdot\rangle_K$。则 $H\times K$ 在通常的向量空间直和运算与内积

$$
\langle(x,\xi),(y,\eta)\rangle_{H\times K}:=\langle x,y\rangle_H+\langle\xi,\eta\rangle_K \tag{3.9}
$$

（$x,y\in H$，$\xi,\eta\in K$）下是内积空间。

此时记 $H\oplus K$（或简记 $H+K$）表示内积空间 $(H\times K,\langle\cdot,\cdot\rangle_{H\times K})$。当 $H$、$K$ 都完备时，称 $H\oplus K$ 为 Hilbert 空间 $H$ 与 $K$ 的**直和**。

**证明** 对 $x,y\in H$、$\xi,\eta\in K$ 逐条验证：

(i) $\langle(x,\xi),(x,\xi)\rangle_{H\times K}=\langle x,x\rangle_H+\langle\xi,\xi\rangle_K\ge0$；若其为零，则 $\langle x,x\rangle_H=0$ 且 $\langle\xi,\xi\rangle_K=0$，故 $x=0$、$\xi=0$，即 $(x,\xi)=0$。

(ii) 对任意 $\tilde x\in H$、$\tilde\xi\in K$，

$$
\begin{aligned}
\langle(x,\xi)+(\tilde x,\tilde\xi),(y,\eta)\rangle_{H\times K}
&=\langle x+\tilde x,y\rangle_H+\langle\xi+\tilde\xi,\eta\rangle_K\\
&=\bigl[\langle x,y\rangle_H+\langle\xi,\eta\rangle_K\bigr]+\bigl[\langle\tilde x,y\rangle_H+\langle\tilde\xi,\eta\rangle_K\bigr]\\
&=\langle(x,\xi),(y,\eta)\rangle_{H\times K}+\langle(\tilde x,\tilde\xi),(y,\eta)\rangle_{H\times K}.
\end{aligned}
$$

(iii) 对任意 $\alpha\in\mathbb K$，

$$
\langle\alpha(x,\xi),(y,\eta)\rangle_{H\times K}=\langle\alpha x,y\rangle_H+\langle\alpha\xi,\eta\rangle_K=\alpha\langle x,y\rangle_H+\alpha\langle\xi,\eta\rangle_K=\alpha\langle(x,\xi),(y,\eta)\rangle_{H\times K}.
$$

(iv) $\langle(x,\xi),(y,\eta)\rangle_{H\times K}=\langle x,y\rangle_H+\langle\xi,\eta\rangle_K=\langle y,x\rangle_H+\langle\eta,\xi\rangle_K=\langle(y,\eta),(x,\xi)\rangle_{H\times K}$。

故 (3.9) 是内积。由 (i)，诱导范数满足

$$
\lVert(x,\xi)\rVert_{H\times K}=\sqrt{\lVert x\rVert_H^2+\lVert\xi\rVert_K^2}.
$$

当 $H$、$K$ 都是 Hilbert 空间时，$H\times K$ 在此范数下显然完备。$\square$

> **要点** 直和范数不是 $\max$ 也不是 $\lVert\cdot\rVert_H+\lVert\cdot\rVert_K$，而是**欧氏合成** $\sqrt{\lVert x\rVert_H^2+\lVert\xi\rVert_K^2}$。这一选择由 (3.9) 的「内积相加」唯一决定 —— 换句话说，**是先要内积，范数才被逼成这个形式**。这与 Banach 空间直和（三种范数等价、任选其一）形成对照，也再次说明内积结构比范数结构更「刚性」。

### 3.1.9 练习

(0) 证明由内积定义的 $\lVert x\rVert=\sqrt{\langle x,x\rangle}$ 是范数。
(1) 证明 $(C[0,1],\lVert\cdot\rVert_\infty)$ 不是内积空间。
(2) $L^p[0,1]$（$1\le p<\infty$，$p\ne2$）不是 Hilbert 空间。
(3) 设 $AC_0$ 为一切满足 $x(0)=0$ 且 $x'\in L^2[0,1]$ 的绝对连续函数 $x:[0,1]\to\mathbb K$ 之集。定义

   $$
   \langle x,y\rangle=\int_0^1x'(t)\overline{y'(t)}\,\mathrm dt,\qquad\forall x,y\in AC_0 .
   $$

   证明 $AC_0$ 是 Hilbert 空间。
(4) 设 $a=(a_n)_{n\in\mathbb N}$ 是正数列。定义

   $$
   \ell^2_a=\Bigl\{(x_n)_{n\in\mathbb N}:\ x_n\in\mathbb C,\ \sum_{n=1}^\infty a_n|x_n|^2<\infty\Bigr\},
   $$

   并在 $\ell^2_a$ 上定义内积 $\langle x,y\rangle=\sum_{n=1}^\infty a_nx_n\overline{y_n}$（$x=(x_n),y=(y_n)\in\ell^2_a$）。证明 $\ell^2_a$ 是 Hilbert 空间。
(5) 实 $\ell^2$ 中一切满足 $|x_n|\le\dfrac1n$（$\forall n\in\mathbb N$）的向量 $x=(x_n)_{n\in\mathbb N}$ 之集称为 **Hilbert 立方体**。证明它在 $\ell^2$ 中紧。
(6) 设 $H$ 是内积空间。
   (a) 若 $a,b\in H\setminus\{0\}$，$a'=\dfrac{a}{\lVert a\rVert^2}$、$b'=\dfrac{b}{\lVert b\rVert^2}$，则

   $$
   \lVert a'-b'\rVert=\frac{\lVert a-b\rVert}{\lVert a\rVert\lVert b\rVert}.
   $$

   (b) 证明 **Ptolemy 不等式**：对一切 $a,b,c,d\in H$，

   $$
   \lVert a-c\rVert\lVert b-d\rVert\le\lVert a-b\rVert\lVert c-d\rVert+\lVert b-c\rVert\lVert a-d\rVert.
   $$

(7) 若 $\lVert\cdot\rVert$ 是向量空间 $V$ 上由内积诱导的范数，证明它满足 **Apollonius 不等式**：对一切 $x,y,z\in V$，

   $$
   \lVert z-x\rVert^2+\lVert z-y\rVert^2=\frac12\lVert x-y\rVert^2+2\Bigl\lVert z-\frac12(x+y)\Bigr\rVert^2.
   $$

(8) 若 $H$ 是 Hilbert 空间、$M$ 是 $H$ 的闭子空间，证明商空间 $H/M$ 也是 Hilbert 空间。

## 练习答案

**(0)** 正定性：$\langle x,x\rangle\ge0$，且 $=0\iff x=0$。齐性：$\langle\lambda x,\lambda x\rangle=\lvert\lambda\rvert^2\langle x,x\rangle$。三角不等式：

$$
\lVert x+y\rVert^2=\lVert x\rVert^2+2\operatorname{Re}\langle x,y\rangle+\lVert y\rVert^2\le\lVert x\rVert^2+2\lVert x\rVert\lVert y\rVert+\lVert y\rVert^2=\bigl(\lVert x\rVert+\lVert y\rVert\bigr)^2,
$$

中间一步用 Cauchy–Schwarz。$\square$

**(1)** 取 $x(t)\equiv1$、$y(t)=t$。则 $\lVert x\rVert_\infty=\lVert y\rVert_\infty=1$、$\lVert x+y\rVert_\infty=2$、$\lVert x-y\rVert_\infty=1$，故

$$
\lVert x+y\rVert_\infty^2+\lVert x-y\rVert_\infty^2=4+1=5\ne2(1+1)=2\bigl(\lVert x\rVert_\infty^2+\lVert y\rVert_\infty^2\bigr).
$$

平行四边形律不成立，故 $\lVert\cdot\rVert_\infty$ 不由内积诱导（定理 3.1.10 的逆否）。$\square$

**(2)** 取 $x=\mathbf 1_{[0,1/2]}$、$y=\mathbf 1_{[1/2,1]}$。则 $\lVert x\rVert_p=\lVert y\rVert_p=2^{-1/p}$，而 $x+y=\mathbf 1_{[0,1]}$、$x-y=\pm\mathbf 1_{[0,1]}$，故 $\lVert x\pm y\rVert_p=1$。平行四边形律要求

$$
1+1=2\bigl(2^{-2/p}+2^{-2/p}\bigr)=2^{2-2/p},
$$

即 $2=2^{2-2/p}$，等价于 $p=2$。故 $p\ne2$ 时 $L^p[0,1]$ 不是内积空间，从而不是 Hilbert 空间。$\square$

> **注** 这一题不能只对**一个** $p$ 举反例就完事 —— 要得到一个与 $p$ 无关的判据，必须让等式化为 $p=2$ 这一个条件。

**(3)** 令 $D:AC_0\to L^2[0,1]$、$Dx=x'$。**线性**显然。**满射**：对 $g\in L^2[0,1]$，令 $x(t)=\int_0^tg$，则 $x$ 绝对连续、$x(0)=0$、$x'=g$ a.e.，故 $Dx=g$。**单射**：$x'=0$ a.e. 且 $x$ 绝对连续、$x(0)=0$ $\Rightarrow x\equiv0$。**等距**：

$$
\lVert Dx\rVert_{L^2}^2=\int_0^1\lvert x'\rvert^2=\langle x,x\rangle=\lVert x\rVert^2 .
$$

故 $AC_0$ 与 $L^2[0,1]$ 等距同构，而后者完备，故 $AC_0$ 完备；内积由 $L^2$ 搬过来，故是 Hilbert 空间。$\square$

**(4)** 令 $U:\ell^2_a\to\ell^2$、$U(x)=(\sqrt{a_n}x_n)_{n\ge1}$。$a_n>0$ 保证 $U$ 是线性双射（逆为 $(y_n)\mapsto(y_n/\sqrt{a_n})$），且

$$
\lVert Ux\rVert_{\ell^2}^2=\sum_na_n\lvert x_n\rvert^2=\lVert x\rVert^2,\qquad \langle Ux,Uy\rangle=\sum_na_nx_n\overline{y_n}=\langle x,y\rangle .
$$

故 $U$ 是等距同构，$\ell^2_a$ 完备且内积由 $\ell^2$ 搬来，是 Hilbert 空间。$\square$

**(5)** 记 $C=\{x\in\ell^2:\lvert x_n\rvert\le1/n\}$。

**闭。** 每个坐标泛函 $x\mapsto x_n$ 连续，故 $C=\bigcap_n\{x:\lvert x_n\rvert\le1/n\}$ 是闭集之交。

**全有界。** 给定 $\varepsilon>0$，取 $N$ 使 $\sum_{n>N}n^{-2}<\varepsilon^2/4$。投影 $\pi_N(C)\subseteq\prod_{n\le N}[-1/n,1/n]$ 是 $\mathbb R^N$ 中的有界闭集，故紧，取有限 $\varepsilon/2$-网 $y^{(1)},\dots,y^{(m)}$；把它们补零看成 $\ell^2$ 中的点。对任意 $x\in C$，取 $j$ 使 $\sum_{n\le N}\lvert x_n-y^{(j)}_n\rvert^2<\varepsilon^2/4$，而 $y^{(j)}$ 的尾巴为零，故

$$
\lVert x-y^{(j)}\rVert^2\le\frac{\varepsilon^2}{4}+\sum_{n>N}\lvert x_n\rvert^2\le\frac{\varepsilon^2}{4}+\sum_{n>N}\frac{1}{n^2}<\frac{\varepsilon^2}{2}<\varepsilon^2 .
$$

故 $C$ 全有界。$\ell^2$ 完备、$C$ 闭故完备，完备 + 全有界 $\Rightarrow$ 紧。$\square$

> **要点** 「有界闭不紧」是无穷维的通病（定理 2.4.8），$C$ 之所以紧，是因为它的各个坐标被 $1/n$ 压得**越来越小**，尾巴整体可弃。这就是「紧性靠尾巴衰减，不靠有界」。

**(6)(a)** 直接展开：

$$
\Bigl\lVert\frac{a}{\lVert a\rVert^2}-\frac{b}{\lVert b\rVert^2}\Bigr\rVert^2=\frac{1}{\lVert a\rVert^2}+\frac{1}{\lVert b\rVert^2}-\frac{2\operatorname{Re}\langle a,b\rangle}{\lVert a\rVert^2\lVert b\rVert^2}=\frac{\lVert a\rVert^2+\lVert b\rVert^2-2\operatorname{Re}\langle a,b\rangle}{\lVert a\rVert^2\lVert b\rVert^2}=\frac{\lVert a-b\rVert^2}{\lVert a\rVert^2\lVert b\rVert^2}.
$$

开方即得。$\square$

**(6)(b)** 用 (a) 加「反演」。固定 $a$，把点 $x\ne a$ 映成

$$
x^*=a+\frac{x-a}{\lVert x-a\rVert^2},
$$

则对 $x,y\ne a$，由 (a)（取 $u=x-a$、$v=y-a$）

$$
\lVert x^*-y^*\rVert=\frac{\lVert x-y\rVert}{\lVert x-a\rVert\lVert y-a\rVert}.
$$

对三点 $b^*,c^*,d^*$ 用三角不等式：

$$
\lVert b^*-d^*\rVert\le\lVert b^*-c^*\rVert+\lVert c^*-d^*\rVert ,
$$

代入上式并乘以 $\lVert b-a\rVert\lVert c-a\rVert\lVert d-a\rVert$：

$$
\lVert b-d\rVert\lVert c-a\rVert\le\lVert b-c\rVert\lVert d-a\rVert+\lVert c-d\rVert\lVert b-a\rVert ,
$$

这正是所求（左端即 $\lVert a-c\rVert\lVert b-d\rVert$，右端两项即 $\lVert b-c\rVert\lVert a-d\rVert+\lVert c-d\rVert\lVert a-b\rVert$）。若某个点等于 $a$，上式两边可直接逐项核对，结论同样成立。$\square$

**(7)** 在平行四边形律 $\lVert u+v\rVert^2+\lVert u-v\rVert^2=2\lVert u\rVert^2+2\lVert v\rVert^2$ 中取 $u=z-x$、$v=z-y$：

$$
u+v=2z-x-y=2\Bigl(z-\frac{x+y}{2}\Bigr),\qquad u-v=y-x .
$$

代入即得

$$
\lVert z-x\rVert^2+\lVert z-y\rVert^2=\frac12\lVert x-y\rVert^2+2\Bigl\lVert z-\frac12(x+y)\Bigr\rVert^2 .
$$

$\square$

> **要点** Apollonius 等式就是「中线性」的另一种写法：它等于说「到两定点的距离平方和只依赖于中点与距离」，四边形四条边与两条对角线的关系（平行四边形律）由此而来。

**(8)** **完备。** $M$ 闭，$H/M$ 是 Banach 空间（2.4.5 练习 7 的反向命题；也可直接用商范数验证 Cauchy 列）。具体地，若 $([x_n])$ Cauchy，取 $u_n\in M$ 使 $\lVert x_n+u_n\rVert$ 逼近 $\lVert[x_n]\rVert$，则 $(x_n+u_n)$ Cauchy，收敛到 $x$，于是 $[x_n]\to[x]$。

**平行四边形律。** 由 3.3 练习 6，$H/M$ 与 $M^\perp$ 等距同构，而 $M^\perp$ 是 Hilbert 空间（闭子空间），故 $H/M$ 满足平行四边形律；内积可取

$$
\langle[x],[y]\rangle_{H/M}=\langle P^\perp x,P^\perp y\rangle .
$$

故 $H/M$ 是 Hilbert 空间。$\square$

> **要点** 商空间丢掉的是「$M$ 方向的信息」，$M^\perp$ 恰好是那份信息的正交补 —— 等距同构把「商」翻译成「投影」。

---

## 订正

1. **定理 3.1.10 证明中 (3.6) 的系数印错了。** 原文在两步平行四边形律之后写

   $$
   \frac14\bigl(\lVert x+y\rVert^2+\lVert z\rVert^2-\lVert x-y\rVert^2-\lVert z\rVert^2\bigr)=2\operatorname{Re}\langle x,y\rangle,
   $$

   但按 (3.4)，左端整体等于 $\operatorname{Re}\langle x,y\rangle$，不是 $2\operatorname{Re}\langle x,y\rangle$。正确系数是 $\dfrac12$：由两次平行四边形律得到的是
   $\dfrac14\cdot2\bigl(\lVert x+y\rVert^2+\lVert z\rVert^2-\lVert x-y\rVert^2-\lVert z\rVert^2\bigr)$。本文正文已改正。
2. **例 3.1.8 的伴随矩阵定义漏了共轭。** 原文写 $B^*=(b_{ji})$，按定义应为 $B^*=(\overline{b_{ji}})$（共轭转置）。否则 $\operatorname{tr}(B^*A)$ 在复情形不满足内积公理 (1)。本文已改。
3. **例 3.1.8 中 Frobenius 范数的括号错位。** 原文印作 $\lVert(a_{ij}\rVert$，应为 $\lVert(a_{ij})\rVert$。
4. **定义 3.1.4 把三角不等式说成「易验证」。** 该步实为 Cauchy–Schwarz 的直接推论，本文在 3.1.3 补充中补出完整一步。
5. **定理 3.1.10 证明末尾「$\langle x,y\rangle=\overline{\langle y,x\rangle}$ is trivial」不算显然。** 由于此处的内积是**构造**出来的，共轭对称性仍须验证（实部由 (3.4) 对称，虚部由 $\operatorname{Re}\langle iy,x\rangle=-\operatorname{Re}\langle ix,y\rangle$）。本文已补。
6. **记号提示（非排印错误）**：原 PDF 的文字层不携带共轭上划线，故定义 3.1.1(2)、注 3.1.2(3')、例 3.1.5–3.1.8 中凡应在第二变量上的共轭在文本抽取时全部丢失。本文已按本书约定（第一变量线性、第二变量共轭线性）补齐。

---

## 下一步

本章后半段顺着「夹角」这条线继续：**3.2 最佳逼近**（$x$ 到闭子空间的最短距离是否取到）、**3.3 正交分解**（$H=M\oplus M^\perp$）、**3.4 正交基**（把 Schauder 基换成更强的正交基，并回答「Hilbert 空间一定有基」）。

---

系列：[[01-预备·Zorn引理与向量空间]] · [[02-度量空间与拓扑空间]] · [[03-拓扑向量空间]] · [[04-赋范空间]] · [[05-Banach空间]] · [[06-可分Banach空间]] · [[07-Banach空间·完备化与紧性]] · [[09-Hilbert空间·最佳逼近]] · [[10-Hilbert空间·正交分解]] · [[11-Hilbert空间·正交基]] · [[12-连续线性算子·有界性与算子范数]] · [[13-连续线性算子·连续线性泛函与弱拓扑]] · [[14-连续线性算子·对偶空间的例子]] · [[15-连续线性算子·伴随算子]] · [[16-连续线性算子·投影算子]]
