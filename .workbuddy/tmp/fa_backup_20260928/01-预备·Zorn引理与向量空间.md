# 01 预备 · Zorn 引理与向量空间

> 系列：[[02-度量空间与拓扑空间]] · [[03-拓扑向量空间]] · [[04-赋范空间]] · [[05-Banach空间]]
>
> 全书固定记号：$\mathbb K$ 表示实数域 $\mathbb R$ 或复数域 $\mathbb C$。

## 1.1 Zorn 引理

### 1.1.1 偏序与定向集

**定义 1.1.1（偏序、定向偏序）** 设 $I\ne\varnothing$。

- $I\times I$ 的任一子集 $\preceq$ 称为 $I$ 上的**二元关系**；$(x,y)\in\preceq$ 记作 $x\preceq y$。
- 二元关系 $\preceq$ 称为**偏序**，若对一切 $x,y,z\in I$：
  1. （自反）$x\preceq x$；
  2. （反对称）$x\preceq y$ 且 $y\preceq x$ $\Longrightarrow$ $x=y$；
  3. （传递）$x\preceq y$ 且 $y\preceq z$ $\Longrightarrow$ $x\preceq z$。
- 偏序 $\preceq$ 称为**定向的**，若对任意 $x,y\in I$，存在 $z\in I$ 使
  $$x\preceq z\quad\text{且}\quad y\preceq z.$$

**例 1.1.2** 非空集 $X$ 的幂集 $\mathcal P(X)$ 按包含关系 $\subseteq$ 构成定向偏序集：$A\preceq B\iff A\subseteq B$；任取 $A,B$，$A\cup B$ 是公共上界。

**例 1.1.3（字典序）** 在 $\mathbb R^2$ 上定义
$$(x_1,y_1)\preceq(x_2,y_2)\iff x_1<x_2,\ \text{或}\ \bigl(x_1=x_2\ \text{且}\ y_1\le y_2\bigr).$$
这是 $\mathbb R^2$ 上的定向偏序，称为**字典序**。

**例 1.1.4（延拓序）** 设 $X,Y$ 是集合，$\mathcal F$ 为一切满足 $D(f)\subseteq X$、$R(f)\subseteq Y$ 的函数 $f$ 之集。定义
$$f\preceq g\iff D(f)\subseteq D(g)\ \text{且}\ f(x)=g(x),\quad \forall x\in D(f),$$
即 $g$ 是 $f$ 的**延拓**。$(\mathcal F,\preceq)$ 是定向偏序集。

**例 1.1.5** 设 $f$ 是 $[0,1]$ 上的实值函数，$\mathcal A$ 为 $[0,1]$ 上被 $f$ 上方控制的一切阶梯函数 $\varphi$（即 $\varphi(x)\le f(x),\ \forall x\in[0,1]$）之集，按逐点序排列。若 $\varphi,\psi$ 均为阶梯函数，则
$$\varphi\vee\psi:\ x\mapsto\max\{\varphi(x),\psi(x)\}$$
也是阶梯函数，故 $(\mathcal A,\le)$ 是定向集。

### 1.1.2 上界、上确界、极大元、链

**定义 1.1.6** 设 $(X,\preceq)$ 是偏序集，$\varnothing\ne Y\subseteq X$。

- $b\in X$ 称为 $Y$（在 $X$ 中）的**上界**，若对一切 $y\in Y$ 有 $y\preceq b$；此时记 $Y\preceq b$。
- $Y$ 的上界 $z$ 称为 $Y$ 的**上确界**，若对 $Y$ 的任一上界 $b$ 都有 $z\preceq b$。上确界若存在则唯一，记作 $\sup Y$。
- $z\in X$ 称为**极大元**，若
  $$z\preceq a\ \Longrightarrow\ z=a,\qquad \forall a\in X.$$
- $X$ 的**链**是指全序子集 $C\subseteq X$，即 $C$ 中任两元素 $a,b$ 可比较：$a\preceq b$ 或 $b\preceq a$。

**例 1.1.7** 取 $\mathcal P(X)$ 按包含偏序，$\{A_n\}$ 递增（$A_n\subseteq A_{n+1}$）。则 $\mathcal A=\{A_n:n\in\mathbb N\}$ 是 $\mathcal P(X)$ 中的链，且 $\bigcup_{n=1}^\infty A_n$ 是 $\mathcal A$ 的上确界。但 $\bigcup_n A_n$ 未必属于 $\mathcal A$：若属于，则 $\mathcal A$ 有极大元 $\bigcup_n A_n$；若不属于，则 $\mathcal A$ 根本没有极大元。
> 这个例子说明：**链有上界**与**链的极大元存在**是两件不同的事，Zorn 引理正是用前者换后者。

### 1.1.3 两条等价公理

**定理 1.1.8（Zorn 引理）** 设 $(X,\preceq)$ 是偏序集，若 $X$ 的每条非空链 $C\subseteq X$ 都有上界（在 $X$ 中），则 $X$ 有极大元。

**定理 1.1.9（选择公理 AC）** 设 $I,X$ 是两个非空集，对每个 $i\in I$ 给定非空子集 $X_i\subseteq X$。则存在映射 $g:I\to X$，使 $g(i)\in X_i$ 对一切 $i\in I$ 成立。

**定理 1.1.10** 选择公理与 Zorn 引理等价。

> **要点** 定理 1.1.8 与 1.1.9 都是 ZFC 集合论系统的**公理**：在 ZFC 其余公理下既不能证明也不能否定它们。但它们彼此等价（证明见教材附录 A）。

---

## 1.2 向量空间

### 1.2.1 定义与例

设 $\mathbb K$ 为 $\mathbb R$ 或 $\mathbb C$。$\mathbb K$ 上的**向量空间**是集合 $X$，配备向量加法 $(a,b)\mapsto a+b:\ X\times X\to X$ 与数乘 $(\lambda,x)\mapsto\lambda x:\ \mathbb K\times X\to X$，满足：

- **(a)** $(X,+)$ 是交换加法群，即
  - (i) $x+y=y+x$，$x+(y+z)=(x+y)+z$；
  - (ii) 存在唯一的零向量 $0$ 使 $x+0=x$；
  - (iii) 对每个 $x$ 存在唯一的加法逆 $-x$ 使 $x+(-x)=0$。
- **(b)** 对一切 $\lambda,\mu\in\mathbb K$、$x,y\in X$：
  $$1x=x,\qquad \lambda(\mu x)=(\lambda\mu)x,\qquad \lambda(x+y)=\lambda x+\lambda y,\qquad (\lambda+\mu)x=\lambda x+\mu x.$$

**例 1.2.1** 基本例子：

| 空间 | 说明 |
|---|---|
| $\mathbb K^n$ | $n$ 元数组，逐坐标加法与数乘 |
| $\mathbb K^{\mathbb N}$ | $\{(x_1,x_2,\dots):x_i\in\mathbb K\}$，无限维，逐坐标运算 |
| $C[0,1]$ | $[0,1]$ 上一切 $\mathbb K$ 值连续函数，逐点运算 |
| $M_n(\mathbb K)$ | $n\times n$ 矩阵，$(a_{ij})+(b_{ij})=(a_{ij}+b_{ij})$，$\lambda(a_{ij})=(\lambda a_{ij})$ |

### 1.2.2 子空间、线性无关、Hamel 基

对 $A,B\subseteq X$、$\lambda\in\mathbb K$ 记
$$A+B=\{x+y:x\in A,\ y\in B\},\qquad \lambda A=\{\lambda x:x\in A\}.$$

$Y\subseteq X$ 称为 $X$ 的**（线性）子空间**，若 $Y$ 在继承的运算下自身是向量空间。等价条件为
$$Y-Y=Y+(-Y)\subseteq Y,\qquad \lambda Y+\mu Y\subseteq Y\quad(\forall\lambda,\mu\in\mathbb K).$$

**例 1.2.2** $\displaystyle\Bigl\{f\in C[0,1]:\int_0^1 f(t)\,dt=0\Bigr\}$ 是 $C[0,1]$ 的线性子空间。

**定义** 子集 $A\subseteq X$ 称为**线性无关**，若对 $A$ 中任意有限多个互异元素 $v_1,\dots,v_n$，
$$\lambda_1v_1+\dots+\lambda_nv_n=0\ \Longrightarrow\ \lambda_1=\dots=\lambda_n=0.$$
$A$ 的**（线性）张成**是
$$\operatorname{span}A=\Bigl\{\sum_{i=1}^n\lambda_iv_i:\ \lambda_i\in\mathbb K,\ v_i\in A,\ n\in\mathbb N\Bigr\}.$$
子集 $B\subseteq X$ 称为 $X$ 的 **Hamel 基**，若 $B$ 线性无关且张成 $X$。此时每个 $x\in X$ 可**唯一**地写成有限线性组合 $x=\lambda_1b_1+\dots+\lambda_kb_k$。

**定理 1.2.3** 每个非零向量空间 $X$ 都有 Hamel 基。

> **证明思路（见附录定理 A.4.1）** 用 Zorn 引理：取一切线性无关子集之族，按包含偏序（见例 1.1.2、1.1.4 的定向上界构造），极大元即为 Hamel 基。这是 Zorn 引理在泛函分析中最典型的第一次应用。

### 1.2.3 练习

1. 设 $M$ 是 $X$ 的子空间，在 $X$ 上定义等价关系 $x\sim y\iff x-y\in M$。商空间 $X/M$ 是等价类 $[x]=\{x+m:m\in M\}$ 之集。验证按 $[x]+[y]=[x+y]$、$\lambda[x]=[\lambda x]$ 它成为 $\mathbb K$ 上的向量空间，且商映射 $Q:x\mapsto[x]$ 线性。
2. 设 $Z$ 是 $X$ 的线性无关子集，证明 $X$ 有包含 $Z$ 的 Hamel 基。
3. 在例 1.2.2 中，若 $\mathbb K=\mathbb R$，证明每个满足 $\int_0^1 f=0$ 的连续实值函数 $f$ 必在 $[0,1]$ 中某点取零值。

---

## 订正

教材原文有两处排版讹误，本文已按数学含义修正，使用时可对照：

1. §1.2 数乘公理最后一条原文印作 $(\lambda+\mu)x=\lambda x+\mu y$，应为 $(\lambda+\mu)x=\lambda x+\mu x$。
2. §1.1.1 传递性原文印作「$x\preceq y$ 且 $y\preceq z$ $\Rightarrow$ $x\le z$」，应为 $x\preceq z$。

---

系列：[[02-度量空间与拓扑空间]] · [[03-拓扑向量空间]] · [[04-赋范空间]] · [[05-Banach空间]]
