# 09 乘积测度与 Fubini 定理

> 系列：[[01-测度与概率空间]] · [[02-测度的扩张]] · [[03-欧氏空间上的测度与可测函数]] · [[04-概率分布]] · [[05-期望与矩（直观）]] · [[06-可测函数的积分]] · [[07-积分的性质与收敛定理]] · [[08-期望（测度论）]] · [[09-乘积测度与Fubini定理]] · [[10-随机向量与联合分布]] · [[11-独立性]]

§1.7 办一件事：在乘积空间上造出测度，并给出**积分换序**的条件。这两件事分别是「联合分布存在」与「独立性可计算」的技术基础 —— 后面 [[10-随机向量与联合分布]] 与 [[11-独立性]] 全靠它。

本节的重点不在定理陈述，而在**三个条件各自不能少的反例** —— 教材专门给了三个 Examples，考试常拿来判断题。

## 9.1 乘积测度

设 $(X,\mathcal A,\mu_1)$、$(Y,\mathcal B,\mu_2)$，其中 $\mu_1,\mu_2$ 皆为 **$\sigma$-有限测度**。

$$
X\times Y=\big\{(x,y):\ x\in X\ \text{且}\ y\in Y\big\},
$$

$$
\mathcal A\times\mathcal B=\sigma(\mathcal S),\qquad \mathcal S=\big\{A\times B:\ A\in\mathcal A,\ B\in\mathcal B\big\}\quad(\text{矩体集}).
$$

**定理 1.7.1（乘积测度的存在唯一性）**　在 $(X\times Y,\ \mathcal A\times\mathcal B)$ 上**存在唯一**的 $\sigma$-有限测度 $\mu$，使得

$$
\mu(A\times B)=\mu_1(A)\mu_2(B),\qquad \forall A\in\mathcal A,\ B\in\mathcal B . \tag{9.1}
$$

称 $\mu$ 为 $\mu_1,\mu_2$ 的**乘积测度**，记为 $\mu_1\times\mu_2$；称 $(X\times Y,\mathcal A\times\mathcal B,\mu_1\times\mu_2)$ 为 $(X,\mathcal A,\mu_1)$ 与 $(Y,\mathcal B,\mu_2)$ 的**乘积测度空间**。

**证明**　按 [[02-测度的扩张]] 的框架（半代数 $\to$ 代数 $\to$ Carathéodory 扩张），只需证明由 (9.1) 定义的集合函数 $\mu(\cdot)$ 在 $\mathcal S$ 上满足**可数可加性**，即

$$
A\times B=\bigsqcup_{i=1}^{\infty}A_i\times B_i\ \Longrightarrow\ \mu(A\times B)=\sum_{i=1}^{\infty}\mu(A_i\times B_i). \tag{9.2}
$$

先验证 $\mathcal S$ 是半代数：$(A\times B)\cap(A'\times B')=(A\cap A')\times(B\cap B')$ 仍属 $\mathcal S$；补集

$$
(A\times B)^{c}=\big((X\setminus A)\times B\big)\sqcup\big(X\times(Y\setminus B)\big)
$$

是 $\mathcal S$ 中两个不交集的并。所以扩张路线可用。

**对固定的 $x\in A$**，看纤维：

$$
\{x\}\times B=\bigcup_{i:\ x\in A_i}\{x\}\times B_i .
$$

由于 $\{A_i\times B_i\}$ 两两不交，对固定的 $x$，满足 $x\in A_i$ 的那些 $B_i$ 两两不交。令

$$
I(x)=\{i:\ x\in A_i\},
$$

则

$$
\mu_2(B)=\mu_2\Big(\bigsqcup_{i\in I(x)}B_i\Big)=\sum_{i\in I(x)}\mu_2(B_i)=\sum_{i=1}^{\infty}\mathbf 1_{A_i}(x)\,\mu_2(B_i),
$$

最后一个等号因为 $i\notin I(x)$ 时 $\mathbf 1_{A_i}(x)=0$。乘以 $\mathbf 1_A(x)$（$x\in A$ 时它等于 $1$）：

$$
\mathbf 1_A(x)\,\mu_2(B)=\sum_{i=1}^{\infty}\mathbf 1_{A_i}(x)\,\mu_2(B_i),\qquad \forall x\in X. \tag{9.3}
$$

**对 (9.3) 两边关于 $\mu_1$ 积分。** 左边：

$$
\int_X\mathbf 1_A(x)\mu_2(B)\,\mu_1(dx)=\mu_1(A)\mu_2(B)=\mu(A\times B).
$$

右边：各项非负，由单调收敛定理（[[07-积分的性质与收敛定理]] 定理 1.5.7）可逐项积分：

$$
\int_X\Big(\sum_{i=1}^{\infty}\mathbf 1_{A_i}(x)\mu_2(B_i)\Big)\mu_1(dx)=\sum_{i=1}^{\infty}\mu_2(B_i)\int_X\mathbf 1_{A_i}(x)\,\mu_1(dx)=\sum_{i=1}^{\infty}\mu_1(A_i)\mu_2(B_i).
$$

即得 (9.2)。$\square$

> **要点**　这个证明的骨架只有一句话：**先把二重的问题降成一维的恒等式 (9.3)，再对另一维积分。** 而 (9.3) 之所以成立，本质上是「纤维上的不交并」—— 把矩体的不交并限制到固定 $x$ 的垂直线上看。这个「先固定一维、再积分」的手法，就是后面 Tonelli/Fubini 定理的雏形；区别在于这里只对示性函数做了一遍。

**注（有限维乘积）**　由定理 1.7.1 与归纳法，若 $(\Omega_i,\mathcal F_i,\mu_i)$（$i=1,\dots,n$）皆为 $\sigma$-有限测度空间，$\Omega=\Omega_1\times\cdots\times\Omega_n$，则 $\sigma$ 代数 $\mathcal F=\sigma\{A_1\times\cdots\times A_n:\ A_i\in\mathcal F_i\}$ 上存在**唯一**测度 $\mu$ 使

$$
\mu(A_1\times\cdots\times A_n)=\prod_{m=1}^{n}\mu_m(A_m).
$$

当 $(\Omega_i,\mathcal F_i,\mu_i)=(\mathbb R,\mathcal B,\lambda)$ 对一切 $i$ 时，得到的是 $\mathbb R^{n}$ 的 Borel 子集上的 **Lebesgue 测度**。

这正好呼应 [[03-欧氏空间上的测度与可测函数]]：$\mathbb R^n$ 上的 Lebesgue 测度既可以由「$\mathcal A_4$ 上半开区间的体积函数 $+$ Carathéodory 扩张」造出来，也可以由「一维 Lebesgue 测度的 $n$ 重乘积」造出来，**两条路给出同一个测度**（唯一性保证）。

**注（$n=\infty$）**　无穷维乘积概率空间

$$
\Big(\prod_{i\in\mathbb N}\Omega_i,\ \prod_{i\in\mathbb N}\mathcal F_i,\ \prod_{i\in\mathbb N}\mu_i\Big)
$$

的存在性见教材 **Theorem 2.1.21**。这是随机变量序列、随机过程（即 $\mathbb T=\mathbb N$ 或 $[0,\infty)$ 上的过程）得以存在的地基。

> **补充（为什么无穷维不能照搬上面的证明）**　有限维的证明靠「矩体的不交并限制到一条纤维上仍是不交并」，这对有限个因子成立；但无穷维的柱集（cylinder set）是**有限个坐标上的限制**的逆像，柱集类只构成**半代数而非代数**，且 $\sigma$ 代数不是「柱集的有限并」而是其生成者。Kolmogorov 扩张定理（Th 2.1.21）的处理办法是：先在柱集上由相容的有限维分布定义集合函数，验证其在**代数**上可数可加，再用 Carathéodory 扩张。**关键难点是相容性条件**（对指标集的任一排列、以及对取边缘分布，有限维分布必须自洽），教材在 §2.1 讲。

## 9.2 Fubini 定理

设 $(X\times Y,\mathcal A\times\mathcal B,\mu_1\times\mu_2)$，$f\in\mathcal A\times\mathcal B$（即 $f$ 关于乘积 $\sigma$ 代数可测）。

**定理 1.7.2（Fubini 定理）**　若 $f\ge0$ 或 $\displaystyle\int\lvert f\rvert\,d\mu<\infty$，则

$$
\int_X\int_Y f(x,y)\,\mu_2(dy)\,\mu_1(dx)=\int_{X\times Y}f\,d\mu=\int_Y\int_X f(x,y)\,\mu_1(dx)\,\mu_2(dy). \tag{*}
$$

（$f\ge0$ 的版本常单列为 **Tonelli 定理**；$f$ 可积的版本是 Fubini。）

**条件满足时，另有可测性结论**：

$$
x\longmapsto\int_Y f(x,y)\,\mu_2(dy)\ \text{为}\ \mathcal A\text{-可测};\qquad y\longmapsto\int_X f(x,y)\,\mu_1(dx)\ \text{为}\ \mathcal B\text{-可测}.
$$

这条不是附注而是**换序等式成立的前提**：若内层积分不可测，外层的 $\int_X(\cdot)\mu_1(dx)$ 根本没有定义。所以严格地说它应该并进定理正文。

> **补充（证明骨架）**　与积分的定义同构：
> 1. $f=\mathbf 1_{A\times B}$：两边都等于 $\mu_1(A)\mu_2(B)$，即 (9.1)；
> 2. $f=\mathbf 1_E$（$E\in\mathcal A\times\mathcal B$）：令 $\mathcal L=\{E:\ \text{换序等式对}\ \mathbf 1_E\ \text{成立}\}$，由 1 知 $\mathcal S\subseteq\mathcal L$，再验证 $\mathcal L$ 是 $\lambda$ 系（对真差与递增并封闭），由 $\pi$-$\lambda$ 定理（见 [[11-独立性]]）得 $\sigma(\mathcal S)=\mathcal A\times\mathcal B\subseteq\mathcal L$；
> 3. 非负简单函数、非负可测函数：线性性与单调收敛；
> 4. 一般可积函数：拆正负部。
>
> 第 2 步的 $\pi$-$\lambda$ 论证是整条链上唯一不平凡的地方 —— 它把「矩体上成立」提升到「整个乘积 $\sigma$ 代数上成立」。

## 9.3 三个条件都不能少

Fubini 定理的条件有三条：$f\ge0$ 或 $\int\lvert f\rvert\,d\mu<\infty$；$\mu_2$（及 $\mu_1$）$\sigma$-有限；$f$ 关于乘积 $\sigma$ 代数可测。教材给了三个反例，**一条件一反例**。

### 1° 不能没有「$f\ge0$ 或 $\int\lvert f\rvert\,d\mu<\infty$」

**例 1.7.5**　取 $X=Y=\{1,2,\dots\}$，$\mathcal A=\mathcal B=$ 全体子集，$\mu_1=\mu_2=$ **计数测度**。对 $m\ge1$，令

$$
f(m,m)=1,\qquad f(m+1,m)=-1,\qquad\text{其它 }f(m,n)=0 .
$$

**断言**

$$
\sum_{m}\sum_{n}f(m,n)=1,\qquad\text{但}\qquad \sum_{n}\sum_{m}f(m,n)=0 .
$$

**验算（外层对 $m$、内层对 $n$）**　固定 $m$：$f(m,n)$ 非零只可能在 $n=m$（值 $1$）与 $n=m-1$（值 $-1$，因为 $f((m-1)+1,m-1)=f(m,m-1)=-1$）。

- $m=1$：只有 $n=1$ 一项，和为 $1$；
- $m\ge2$：$n=m$ 给 $1$、$n=m-1$ 给 $-1$，和为 $0$。

故 $\sum_m\sum_n f(m,n)=1+0+0+\cdots=1$。

**验算（外层对 $n$、内层对 $m$）**　固定 $n$：$f(m,n)$ 非零只可能在 $m=n$（值 $1$）与 $m=n+1$（值 $-1$）。两项相消，和为 $0$ 对一切 $n$。故 $\sum_n\sum_m f(m,n)=0$。

**病根**　$\sum_{m,n}\lvert f(m,n)\rvert=\infty$（非零项有无穷多个，每个绝对值 $1$），**不满足可积性**；同时 $f$ 也不非负。所以两个条件都不满足，定理不适用是正当的。

> **要点**　这条反例的结构是「**把质量沿对角线平移**」：$+1$ 在 $(m,m)$、$-1$ 在 $(m+1,m)$。沿行求和时每行只剩一个 $+1$（第一行）或相消（其余行），总和 $1$；沿列求和时每列恰好一对相消，总和 $0$。**读的时候最容易把内外层看反**，建议按上面的「固定一个、枚举另一个」逐项核对。

### 2° 不能没有「$\mu_2$ 是 $\sigma$-有限的」

**例 1.7.7**　取 $X=(0,1)$，$\mathcal A=$ Borel 集，$\mu_1=$ **Lebesgue 测度**；$Y=(0,1)$，$\mathcal B=$ **全体子集**，$\mu_2=$ **计数测度**。令

$$
f(x,y)=\begin{cases}1,&x=y,\\0,&x\ne y.\end{cases}
$$

**验算**　固定 $x$：$\{y:\ f(x,y)=1\}=\{x\}$ 是单点集，计数测度为 $1$，故 $\displaystyle\int_Y f(x,y)\mu_2(dy)=1$ 对一切 $x$，于是

$$
\int_X\int_Y f(x,y)\,\mu_2(dy)\,\mu_1(dx)=\int_X 1\,dx=1 .
$$

固定 $y$：$\{x:\ f(x,y)=1\}=\{y\}$ 是单点集，Lebesgue 测度为 $0$，故 $\displaystyle\int_X f(x,y)\mu_1(dx)=0$ 对一切 $y$，于是

$$
\int_Y\int_X f(x,y)\,\mu_2(dx)\,\mu_2(dy)=\int_Y 0\,d(\text{计数})=0 .
$$

**病根**　$\mu_2=$ 计数测度在 $(0,1)$ 上**不是 $\sigma$-有限的**：$(0,1)$ 不可数，任何可数并的可数集仍是可数集，覆盖不了 $(0,1)$。而 $\mu_1$ 是可数的单点集给 $0$、$\mu_2$ 给 $1$ —— **两种测度对「单点集」的赋值差了一个量级**，这就是换序不成立的原因。

> **要点**　这条反例的病根常被误认为「$f$ 不可测」。其实 $f=\{(x,x):x\in(0,1)\}$ 是**对角线**，它关于 $\mathcal A\times\mathcal B$ 可测（$\mathcal B$ 是全体子集，任何集合都可测）。真正失效的是 $\sigma$-有限性。

### 3° 不能没有「$f$ 关于乘积 $\sigma$ 代数可测」

**例 1.7.8**　讲义只给了编号，没有展开。

**补充（标准的不可测反例）**　取 $X=Y=(0,1)$，$\mu_1=\mu_2=$ Lebesgue 测度，$V\subset(0,1)$ 为 **Vitali 不可测集**（用选择公理构造：在 $\mathbb Q$ 作用的等价类中各取一个代表）。令

$$
f(x,y)=\begin{cases}1,&x-y\in\mathbb Q\ \text{且}\ x\in V+\mathbb Q,\\0,&\text{其它}.\end{cases}
$$

更简单的版本：取 $f=\mathbf 1_{V\times\{y_0\}}$，其中 $V$ 不可测、$y_0$ 固定。则
- 固定 $x$：$f(x,y)=\mathbf 1_{\{y_0\}}(y)\mathbf 1_V(x)$，故 $\int_Y f(x,y)dy=\mathbf 1_V(x)\cdot0=0$（单点的 Lebesgue 测度为 $0$），外层的 $\int_X$ 也存在且为 $0$；
- 固定 $y$：$y\ne y_0$ 时 $\int_X f\,dx=0$；$y=y_0$ 时 $\int_X f\,dx=\int\mathbf 1_V\,dx$，而 $V$ 不可测，**这个积分不存在**。

所以右边那个累次积分根本没有定义。病根是 $V\times\{y_0\}$ 不属于 $\mathcal B\times\mathcal B$（否则它的 $x$-截面 $V$ 必须属于 $\mathcal B$），$f$ 不可测。

讲义的例 1.7.8 只写了标题。**这一条的要点是：可测性不是技术条件，它决定了内层积分是否是一个可测函数、从而外层积分是否有意义。**

## 9.4 推论：乘积测度为零的刻画

**推论**　设 $(X\times Y,\mathcal A\times\mathcal B,\mu_1\times\mu_2)$，$E\in\mathcal A\times\mathcal B$。则

$$
\mu_1\times\mu_2(E)=0
$$

的**充要条件**是

$$
\mu_2(E_x)=0\quad\mu_1\text{-a.e. }x,\qquad\text{且}\qquad \mu_1(E_y)=0\quad\mu_2\text{-a.e. }y,
$$

其中

$$
E_x=\{y:\ (x,y)\in E\}\quad(\forall x\in X),\qquad E_y=\{x:\ (x,y)\in E\}\quad(\forall y\in Y).
$$

**证明**　对 $\mathbf 1_E\ge0$ 用 Tonelli 定理（定理 1.7.2 的 $f\ge0$ 版本）：

$$
\mu_1\times\mu_2(E)=\int_X\mu_2(E_x)\,\mu_1(dx)=\int_Y\mu_1(E_y)\,\mu_2(dy).
$$

两个被积函数都非负，故积分（可能为 $+\infty$）为零当且仅当被积函数 a.e. 为零。$\square$

> **要点**　这条推论是后面判「两个随机变量独立 $\Rightarrow$ 联合分布绝对连续」之类的工具，也是「零测集的截口几乎都是零测集」这句常用话的严格版本。**注意它依赖 $\sigma$-有限性**：例 1.7.7 里 $E=\{(x,x)\}$ 满足 $\mu_2(E_x)=1\ne0$，所以乘积测度不为零（实际上为 $+\infty$），与推论不矛盾 —— 但若把推论用在非 $\sigma$-有限的情形，会出现「每个截口都零测、乘积测度却非零」的假象。

## 订正

1. **「直接应用 Th 1.1.9（半代数）」——这个编号不存在。**　讲义第 67 页写「Th 1.7.1 的证明：直接应用 Th 1.1.9（半代数）」，但 §1.1 里承担「半代数上的测度可扩张到 $\sigma(\mathcal S)$」这一步的定理编号是 **1.1.2**（Stieltjes 测度扩张，见 [[02-测度的扩张]]）或教材的 **Theorem A.1.3**（Carathéodory 扩张定理，讲义第 8 页截图）。§1.1 中出现过的编号只有 1.1.1、1.1.2、1.1.4、1.1.6、1.1.11，**没有 1.1.9**。本笔记按定理 1.1.2 的证明框架（半代数 $\to$ 代数 $\to$ Carathéodory）引用。
2. **定理 1.7.1 的证明讲义只写到「$\mathbf 1_A(x)\mu(B)=\sum\mathbf 1_{A_i}(x)\mu(B_i)$」就停了。**　最后一步「对 $x$ 积分得可数可加性」没有写。已在正文补出，并说明这里用的是**单调收敛定理**（各项非负，可逐项积分），不需要 Fubini —— 否则就成了循环论证（用 Fubini 证乘积测度存在）。
3. **「内层积分的可测性」被排在定理之后当附注。**　讲义第 69 页把它写在右边一栏「Fubini 定理条件满足下」，形式上像补充说明。实际上它是定理成立的前提（否则外层积分无定义），应并入定理陈述。已在正文点明。
4. **例 1.7.7 的笔误。**　讲义写「$\int_Y\int_X f(x,y)\mu_1(dx)\mu_2(dx)=0$」，第二个测度应为 $\mu_2(dy)$（外层对 $y$ 积分）。同一行上方的公式 $\int_X f(x,y)\mu_1(dx)=0$ 是对的。已在正文改正。
5. **例 1.7.8 只有标题、没有内容。**　讲义第 69 页写「3° $f$ 可测性 / 例 1.7.8」后即翻页结束，没有给出反例。已在**补充**中给出标准构造（Vitali 不可测集），并说明病根。
6. **例 1.7.5 的「$f(m,m)=1$、$f(m+1,m)=-1$」与「$\sum_m\sum_n=1$」需要逐项核对才可信。**　讲义只给结论。正文给出了完整的「固定一个、枚举另一个」验算 —— 这条反例的内外层极易看反，**建议照正文的验算表自己复算一遍再记结论**。
7. **注（$n=\infty$）的「Th 2.1.21」在本讲义范围内没有内容。**　讲义只写「见 Th 2.1.21」，教材 §2.1 才讲。本笔记在**补充**里给出了 Kolmogorov 扩张定理的路线（柱集上的相容有限维分布 $\to$ 代数上可数可加 $\to$ Carathéodory 扩张），并指出**相容性条件是难点**。这条是随机过程存在性的地基，等有 §2.1 的材料再补完整证明。
