# 07 乘积测度与 Fubini 定理

> 系列：[[01-概率空间]] · [[02-随机变量与分布]] · [[03-随机变量与随机过程]] · [[04-可测函数与积分]] · [[05-积分的性质与收敛定理]] · [[06-期望]] · [[07-乘积测度与Fubini定理]]

这一节解决「高维怎么积」。逻辑上只有两步：先造出乘积测度 $\mu_1\times\mu_2$（定理 1.7.1），再证明重积分可以换序（定理 1.7.2，Fubini）。两个结论都靠同一个套路：**半代数上的可数可加性 + Carathéodory 扩张**。关键是记住 Fubini 的三个条件 —— 草稿把这三个条件的**反例**全部列了出来，它们是本章最容易被忽略也最容易考的部分。

## 1.7.1 乘积测度

设 $(X,\mathcal A,\mu_1)$、$(Y,\mathcal B,\mu_2)$，其中 $\mu_1,\mu_2$ 皆为 **$\sigma$-有限测度**。

$$
X\times Y=\big\{(x,y):\ x\in X\ \text{且}\ y\in Y\big\},
$$

$$
\mathcal A\times\mathcal B=\sigma(\mathcal S),\qquad \mathcal S=\big\{A\times B:\ A\in\mathcal A,\ B\in\mathcal B\big\}\quad(\text{矩体集}).
$$

**定理 1.7.1（乘积测度的存在唯一性）**　在 $(X\times Y,\ \mathcal A\times\mathcal B)$ 上**存在唯一**的 $\sigma$-有限测度 $\mu$，使得

$$
\mu(A\times B)=\mu_1(A)\mu_2(B).
$$

称 $\mu$ 为 $\mu_1,\mu_2$ 的**乘积测度**，记为 $\mu_1\times\mu_2$；称 $(X\times Y,\mathcal A\times\mathcal B,\mu_1\times\mu_2)$ 为 $(X,\mathcal A,\mu_1)$ 与 $(Y,\mathcal B,\mu_2)$ 的**乘积测度空间**。

**证明**　由定理 1.1.4 的证明框架（半代数 $\to$ 代数 $\to$ Carathéodory 扩张），只需证明 $\mu(\cdot)$ 在 $\mathcal S$ 上满足**可数可加性**，即

$$
A\times B=\bigsqcup_{i=1}^{\infty}A_i\times B_i\ \Longrightarrow\ \mu(A\times B)=\sum_{i=1}^{\infty}\mu(A_i\times B_i),
$$

其中 $A_i\times B_i$ 两两不交。

**第一步：把二重积分拆成可数和。**　对任意 $x\in A$，截面

$$
\{x\}\times B=\bigcup_{x\in A_i}B_i,
$$

这是一个**不交并**：因为 $A_i\times B_i$ 两两不交，对固定的 $x$，至多有一个 $i$ 使 $x\in A_i$；更准确地说，若 $x\in A_i\cap A_j$ 且 $i\ne j$，则 $B_i\cap B_j=\varnothing$。令

$$
I(x)=\{i:\ x\in A_i\},
$$

则 $\{B_i,\ i\in I(x)\}$ 两两不交。于是对固定 $x\in A$，

$$
\mathbf 1_A(x)\,\mu_2(B)=\mu_2\Big(\bigcup_{i\in I(x)}B_i\Big)=\sum_{i\in I(x)}\mu_2(B_i)=\sum_{i=1}^{\infty}\mathbf 1_{A_i}(x)\,\mu_2(B_i),
$$

其中第二个等号用了 $\mu_2$ 的可数可加性，第三个等号把「$i\in I(x)$」写成「$x\in A_i$」。

**第二步：对 $x$ 关于 $\mu_1$ 积分。**　上式两边对 $x$ 关于 $\mu_1$ 积分（用积分换元与线性性）：

$$
\mu_1(A)\,\mu_2(B)=\int_X\mathbf 1_A(x)\mu_2(B)\,\mu_1(dx)=\sum_{i=1}^{\infty}\int_X\mathbf 1_{A_i}(x)\mu_2(B_i)\,\mu_1(dx)=\sum_{i=1}^{\infty}\mu_1(A_i)\mu_2(B_i).
$$

即 $\mu(A\times B)=\sum_i\mu(A_i\times B_i)$，$\mu$ 在 $\mathcal S$ 上可数可加。$\square$

> **要点**　这个证明的漂亮之处在于：**二重级数的可加性被两次一重可加性代替**。第一次用的是 $\mu_2$ 在 $Y$ 上的可数可加性（对每个固定的 $x$），第二次用的是 $\mu_1$ 在 $X$ 上的可数可加性。整个论证不需要 $\mathcal S$ 是 $\sigma$ 代数，只需要它是半代数 —— 这正是 [[01-概率空间]] 定理 1.1.4 那套框架的复用。

**注 (i)（有限多个因子的乘积）**　由定理 1.7.1 与归纳，若 $(\Omega_i,\mathcal F_i,\mu_i),\ i=1,\dots,n$ 是 $\sigma$-有限测度空间、$\Omega=\Omega_1\times\dots\times\Omega_n$，则在由形如 $A_1\times\dots\times A_n\ (A_i\in\mathcal F_i)$ 的集合生成的 $\sigma$ 代数 $\mathcal F$ 上，存在唯一测度 $\mu$ 使

$$
\mu(A_1\times\dots\times A_n)=\prod_{m=1}^{n}\mu_m(A_m).
$$

当 $(\Omega_i,\mathcal F_i,\mu_i)=(\mathbb R,\mathcal R,\lambda)$ 对一切 $i$ 都成立时，得到的即是 **$n$ 维欧氏空间 $\mathbb R^{n}$ 的 Borel 子集上的 Lebesgue 测度**。

**注 (ii)（无穷多个因子的乘积）**　$n=\infty$ 时是**无穷维乘积概率空间** $\Big(\prod_{i\in\mathbb N}\Omega_i,\ \prod_{i\in\mathbb N}\mathcal F_i,\ \prod_{i\in\mathbb N}\mu_i\Big)$，见教材定理 2.1.21。这就是 [[03-随机变量与随机过程]] 中「随机变量序列」与「随机过程」的存在性所在的空间。

## 1.7.2 Fubini 定理

设 $(X\times Y,\ \mathcal A\times\mathcal B,\ \mu_1\times\mu_2)$，$f\in\mathcal A\times\mathcal B$（即 $f$ 关于乘积 $\sigma$ 代数可测）。

**定理 1.7.2（Fubini 定理）**　若 $f\ge0$ 或 $\displaystyle\int|f|\,d\mu<\infty$，则

$$
\int_{X}\int_{Y}f(x,y)\,\mu_2(dy)\,\mu_1(dx)=\int_{X\times Y}f\,d\mu=\int_{Y}\int_{X}f(x,y)\,\mu_1(dx)\,\mu_2(dy). \tag{*}
$$

即**积分可以换序**。

**当 Fubini 定理条件满足时**，还可以得到两个**可测性**结论：

$$
x\longmapsto\int_Y f(x,y)\,\mu_2(dy)\qquad\text{为 }\mathcal A\text{-可测};
$$

$$
y\longmapsto\int_X f(x,y)\,\mu_1(dx)\qquad\text{为 }\mathcal B\text{-可测}.
$$

> **要点**　注意定理的表述里其实含了**三个断言**：(1) 内层积分作为外层变量的函数是可测的；(2) 内层积分存在（可积或非负）；(3) 两次积分相等。写证明时三个都要交代 —— 只验证数值相等是不完整的。$f\ge0$ 的情形通常单列并称为 **Tonelli 定理**：此时不需要可积性，换序与可测性都成立；然后再对 $f=f^{+}-f^{-}$ 用 Tonelli 得到 Fubini。

## 1.7.3 Fubini 定理的三个条件都不能少

草稿把「注意验证条件」明确列为三条，并各配一个反例。

### 1° 不能没有「$f\ge0$ 或 $\int|f|\,d\mu<\infty$」

**例 1.7.5**　设 $X=Y=\{1,2,\dots\}$，$\mathcal A=\mathcal B=$ 全体子集，$\mu_1=\mu_2=$ **计数测度**。对 $m\ge1$，令

$$
f(m,m)=1,\qquad f(m+1,m)=-1,\qquad \text{其它}\ f(m,n)=0.
$$

则

$$
\sum_{m}\sum_{n}f(m,n)=1,\qquad\text{但}\qquad \sum_{n}\sum_{m}f(m,n)=0.
$$

**补充（验算）**
**内层对 $n$ 求和**：固定 $m$，$f(m,n)\ne0$ 的位置是 $n=m$（值 $1$）与 $n=m-1$（值 $-1$，仅当 $m\ge2$）。故

$$
\sum_n f(m,n)=\begin{cases}1,&m=1\\ 1-1=0,&m\ge2\end{cases}\qquad\Longrightarrow\qquad \sum_m\sum_n f(m,n)=1.
$$

**内层对 $m$ 求和**：固定 $n$，$f(m,n)\ne0$ 的位置是 $m=n$（值 $1$）与 $m=n+1$（值 $-1$）。故 $\sum_m f(m,n)=1-1=0$，从而 $\sum_n\sum_m f(m,n)=0$。

两个结果不同，说明**没有可积性（或非负性）保证时换序是错的**。这里 $\iint|f|\,d(\mu_1\times\mu_2)=\infty$（每个 $m$ 贡献 $2$），条件被违反。

### 2° 不能没有「$\mu_2$ 是 $\sigma$-有限的」

**例 1.7.7**　设 $X=(0,1)$，$\mathcal A=$ Borel 集，$\mu_1=$ Lebesgue 测度；$Y=(0,1)$，$\mathcal B=$ 全体子集，$\mu_2=$ **计数测度**。令

$$
f(x,y)=\begin{cases}1,&x=y\\ 0,&x\ne y\end{cases}
$$

则

$$
\int_Y f(x,y)\,\mu_2(dy)=1\quad\text{对一切 }x\qquad\Longrightarrow\qquad \int_X\int_Y f(x,y)\,\mu_2(dy)\,\mu_1(dx)=1,
$$

$$
\int_X f(x,y)\,\mu_1(dx)=0\quad\text{对一切 }y\qquad\Longrightarrow\qquad \int_Y\int_X f(x,y)\,\mu_1(dx)\,\mu_2(dy)=0.
$$

> **要点**　注意这里的**乘数位置**：是两个**外层积分**不同 —— 一次对 $\mu_1$ 积（得 $1$），一次对 $\mu_2$ 积（得 $0$）。原因是 $\mu_2=$ 计数测度在不可数集 $(0,1)$ 上**不是 $\sigma$-有限的**（任何覆盖 $(0,1)$ 的可数并都不可数，故 $\mu_2\equiv\infty$）。反例与 1° 的区别在于：这里**内层积分对每个固定变量都是有限的、可测的**，问题出在把「至多可数条非零对角线」与「不可数条」交换次序上。

### 3° 不能没有「$f$ 关于乘积 $\sigma$ 代数可测」

**例 1.7.8**　存在 $f$ 关于 $x$、$y$ 分别可测，但**不**关于 $\mathcal A\times\mathcal B$ 可测（教材例 1.7.8）。此时内层积分的可测性结论失效。

> **补充（为什么会有这种情形）**　$\mathcal A\times\mathcal B$ 由「矩形」生成，但矩形类**不**包含「$x=y$」这种对角集。要让每个截面的可测性「拼」成整体的可测性，正是 Tonelli 定理的技术核心（对非负 $f$，先用简单函数逼近，逐层验证截面可测）。所以条件 3° 不是形式主义。

## 1.7.4 推论：乘积测度为零的刻画

设 $(X\times Y,\mathcal A\times\mathcal B,\mu_1\times\mu_2)$，$E\in\mathcal A\times\mathcal B$。则

$$
\mu_1\times\mu_2(E)=0
$$

的**充要条件**是

$$
\mu_2(E_x)=0\quad \mu_1\text{-a.e. }x,\qquad\text{且}\qquad \mu_1(E_y)=0\quad \mu_2\text{-a.e. }y,
$$

其中对 $\forall x\in X$、$E_x:=\{y:\ (x,y)\in E\}$；对 $\forall y\in Y$、$E_y:=\{x:\ (x,y)\in E\}$。

**补充（证明）**　对 $f=\mathbf 1_E\ge0$ 应用定理 1.7.2（Tonelli 情形）：

$$
\mu_1\times\mu_2(E)=\int_{X\times Y}\mathbf 1_E\,d\mu=\int_X\mu_2(E_x)\,\mu_1(dx)=\int_Y\mu_1(E_y)\,\mu_2(dy).
$$

被积函数 $\mu_2(E_x)\ge0$ 可测，故第一个积分 $=0$ $\iff$ $\mu_2(E_x)=0$ 对 $\mu_1$-a.e. $x$。同理第二个积分 $=0$ $\iff$ $\mu_1(E_y)=0$ 对 $\mu_2$-a.e. $y$。两边都等于 $\mu_1\times\mu_2(E)$，故等价。

**用途**：这是处理「独立随机变量的事件是否几乎必然不发生」的标准工具。取 $X,Y$ 为两个概率空间、$E$ 为「两坐标都落在某个坏集」时即得零一律型结论。

## 1.7.5 练习

1. 验证例 1.7.5 中 $\iint|f|\,d(\mu_1\times\mu_2)=\infty$，并说明它违反的是 Fubini 的哪一条条件。
2. 在例 1.7.5 里把计数测度换成 $\mu_1=\mu_2=$ 「除有限集外为 $0$ 的测度」的某个 $\sigma$-有限版本，看结论是否改变。
3. 证明：若 $\mu_1,\mu_2$ 都是概率测度，则 $\mu_1\times\mu_2$ 也是概率测度，且 $\mu_1\times\mu_2(X\times Y)=1$。
4. 由定理 1.7.2 推出：$X$ 与 $Y$ 独立的充要条件是，对一切 Borel 集 $A,B$ 有 $P(X\in A,Y\in B)=P(X\in A)P(Y\in B)$；并由此说明此时 $\mu_{X,Y}=\mu_X\times\mu_Y$。
5. 用 1.7.4 的推论证明：若 $X,Y$ 独立、$P(X\in A)>0$、$P(Y\in B)>0$，则 $E=A\times B$ 满足 $\mu_1\times\mu_2(E)>0$。

---
## 订正

1. **引用「Th 1.1.9」不成立。**　草稿证明定理 1.7.1 时写「直接应用 Th 1.1.9（半代数）」，但 1.1 节里承担「半代数上的测度可扩张到 $\sigma(\mathcal S)$」这一步的定理，草稿中一处编号为 Theorem 1.1.4（教材原文截图）、另一处手写为 Th 1.1.2。本笔记统一按**定理 1.1.4 / Carathéodory 扩张定理（Theorem A.1.3）**引用（见 [[01-概率空间]] 的订正节第 2、3 条）。
2. **例 1.7.5 的两个和容易抄反。**　草稿给出「$\sum_m\sum_n f(m,n)=1$ 但 $\sum_n\sum_m f(m,n)=0$」。经逐项验算无误，但**读的时候极易把内外层看反**。本笔记在 1.7.3 的补充里给出了完整的逐项验算。
3. **Fubini 定理的表述中「内层积分的可测性」被省略。**　草稿只写了换序等式 $(*)$，并另起一行列了 $x\mapsto\int_Yf\,\mu_2(dy)$ 是 $\mathcal A$-可测。本笔记把可测性并入定理正文，作为定理的一部分，因为它是换序等式**成立的前提**而不是附注。
4. **「$n=\infty$ 的乘积概率空间」只有一句指向。**　草稿写「见 Th 2.1.21」而未写内容。本笔记保留该指向，并说明它就是随机变量序列/随机过程的存在性空间（[[03-随机变量与随机过程]] 1.3.7）；**该定理的具体内容不在本 PDF 范围内**，引用时需查教材第 2 章。

---
系列：[[01-概率空间]] · [[02-随机变量与分布]] · [[03-随机变量与随机过程]] · [[04-可测函数与积分]] · [[05-积分的性质与收敛定理]] · [[06-期望]] · [[07-乘积测度与Fubini定理]]
