# 数据处理不等式、Fano 与信源编码定理

> 系列：[[01-熵、相对熵与互信息]]
> 本篇按课堂手写稿（3 页）整理，编号自编。记法与 [[01-熵、相对熵与互信息]] 一致：随机变量大写、分布小写、散度双竖线写 $\lVert$。

## 1 马氏链

**定义 1.1（马氏链）** 称随机变量组 $X_1,X_2,\dots,X_n$ 构成马氏链（Markov chain），记作 $X_1\to X_2\to\cdots\to X_n$，若对任意满足 $p(x_1^n)>0$ 的取值

$$
p(x_1^n)=p(x_1)\,p(x_2\mid x_1)\,p(x_3\mid x_2)\cdots p(x_n\mid x_{n-1}),
$$

在不满足 $p(x_1^n)>0$ 的点上该式两端同为 $0$，不必另作规定。等价地，对每个 $i=2,\dots,n$

$$
p(x_i\mid x_1,x_2,\dots,x_{i-1})=p(x_i\mid x_{i-1}),
$$

即「已知现在，过去与将来条件独立」。

> **要点** 定义里那个「其他情况取 $0$」是多余的：若某个 $p(x_i)=0$，则左端 $p(x_1^n)=0$，右端也含零因子。真正要小心的是 $0\cdot\infty$ 型，那是把 $p(x_i\mid x_{i-1})$ 与 $p(x_{i-1})=0$ 同时出现的情形 —— 约定「条件概率在条件事件零概率时不起作用」即可，两端都为 $0$。

**命题 1.2（反向仍为马氏链）** 若 $X_1\to\cdots\to X_n$ 是马氏链，则 $X_n\to X_{n-1}\to\cdots\to X_1$ 也是马氏链。

**证明** 只需验证，在 $p(x_2^n)>0$ 时 $p(x_1\mid x_2,\dots,x_n)=p(x_1\mid x_2)$。由马氏性把联合分布按正向拆开，

$$
p(x_1\mid x_2^n)=\frac{p(x_1^n)}{p(x_2^n)}=\frac{p(x_1)p(x_2\mid x_1)\displaystyle\prod_{i=3}^{n}p(x_i\mid x_{i-1})}{\displaystyle\sum_{x_1'}p(x_1')p(x_2\mid x_1')\prod_{i=3}^{n}p(x_i\mid x_{i-1})}.
$$

分母里的连乘只依赖 $(x_2,\dots,x_n)$，与求和变量 $x_1'$ 无关；分子的连乘同理与 $x_1$ 无关。把它们约掉，得

$$
p(x_1\mid x_2^n)=\frac{p(x_1)p(x_2\mid x_1)}{\sum_{x_1'}p(x_1')p(x_2\mid x_1')}=p(x_1\mid x_2).
$$

$\square$

> **要点** 这个结论是「时间可逆」在条件独立意义下的版本。它保证「下游看上游」与「上游看下游」都有马氏性，于是数据处理不等式可以对称地用在两个方向上（见推论 2.3）。

**命题 1.3（子链与合并）** 设 $X_1\to\cdots\to X_n$ 是马氏链，$A_1,\dots,A_k$ 是 $\{1,\dots,n\}$ 的**按序不交**子集，即 $\max A_j<\min A_{j+1}$。记 $X_{A_j}=(X_i:i\in A_j)$，则

$$
X_{A_1}\to X_{A_2}\to\cdots\to X_{A_k}
$$

也是马氏链。特别地，任取下标 $i_1<i_2<\cdots<i_k$，子序列 $X_{i_1}\to X_{i_2}\to\cdots\to X_{i_k}$ 是马氏链。

**证明** 只用一个手段：把正向乘积按指标集切段，再对不要的变量求和。

记 $m_j=\max A_j$（按序不交保证 $m_1<\cdots<m_k$），并约定 $m_0=0$。把正向分解的因子按右端点分段，第 $j$ 段取

$$
\mathcal F_j=\big\{p(x_i\mid x_{i-1}):m_{j-1}<i\le m_j\big\},
$$

其中 $i=1$ 的因子写成 $p(x_1\mid x_0)=p(x_1)$。每条因子只连接相邻两个下标，故第 $j$ 段涉及的变量全部落在 $[m_{j-1},m_j]$ 内，与更左、更右的段都不相交。

对第 $j$ 段中不属于 $A_j$ 的变量求和（左端点 $x_{m_{j-1}}$ 保留，它是上一段的末位），得

$$
\kappa_j\big(x_{A_j}\mid x_{m_{j-1}}\big)=\sum_{x_i:\,i\in[m_{j-1},m_j]\setminus(A_j\cup\{m_{j-1}\})}\ \prod_{i\in\mathcal F_j}p(x_i\mid x_{i-1}).
$$

由于链式因子逐步归一，$\kappa_j$ 关于 $x_{A_j}$ 是概率核：固定 $x_{m_{j-1}}$ 后对 $x_{A_j}$ 求和得 $1$。$m_k$ 之后的那段尾巴（$i>m_k$）各自对 $x_i$ 求和也都得 $1$。于是

$$
p\big(x_{A_1},\dots,x_{A_j}\big)=\prod_{l\le j}\kappa_l\big(x_{A_l}\mid x_{m_{l-1}}\big),
$$

因为 $l>j$ 的每个 $\kappa_l$ 对 $x_{A_l}$ 求和都为 $1$。两式相除：

$$
p\big(x_{A_j}\mid x_{A_1},\dots,x_{A_{j-1}}\big)=\kappa_j\big(x_{A_j}\mid x_{m_{j-1}}\big)=p\big(x_{A_j}\mid x_{m_{j-1}}\big),
$$

而 $m_{j-1}\in A_{j-1}$，故它只依赖 $x_{A_{j-1}}$。这就是马氏性。取 $A_j=\{i_j\}$ 即得子链情形。$\square$

> **要点** 证明里真正起作用的只有一件事：**链式因子的支撑是局部的**（每条只碰相邻两个下标）。所以「切成按序不交的若干段」之后，段与段之间唯一的接口是端点 $x_{m_{j-1}}$，条件独立就自动成立。这条思路会反复出现 —— AEP 里「典型序列的概率因式分解」、信源编码里「$X^n\to\hat X^n$ 是马氏链」用的都是同一个局部性。

**推论 1.4** 若 $X_1\to\cdots\to X_n$ 是马氏链，$f_1,\dots,f_n$ 是确定性函数，则 $f_1(X_1)\to f_2(X_2)\to\cdots\to f_n(X_n)$ 也是马氏链。

**证明** $f_i(X_i)$ 是 $X_i$ 的函数，给定 $X_i$ 后它没有额外随机性，条件独立性原样传递。$\square$

## 2 数据处理不等式

**定理 2.1（数据处理不等式，DPI）** 若 $X\to Y\to Z$ 是马氏链，则

$$
I(X;Z)\le I(X;Y).
$$

**证明** 由互信息的链式法则，把 $I(X;Y,Z)$ 按两个方向展开：

$$
I(X;Y,Z)=I(X;Y)+I(X;Z\mid Y)=I(X;Y),
$$

第一个等号处用马氏性得 $I(X;Z\mid Y)=0$；另一个方向

$$
I(X;Y,Z)=I(X;Z)+I(X;Y\mid Z)\ge I(X;Z),
$$

末步用 $I(X;Y\mid Z)\ge0$。两式比较即得。$\square$

**补充（一条恒等式）** 证明里两个方向的展开给出的是同一个量 $I(X;Y,Z)$，这不是巧合：它说明「$I(X;Y)$ 与 $I(X;Z)$ 的差正好是 $I(X;Y\mid Z)$」。

$$
I(X;Y)-I(X;Z)=I(X;Y\mid Z)\ge0 .
$$

这条恒等式在等号条件的讨论里反复出现，值得单独记住。

> **要点（怎么读这条不等式）** 「处理数据不能增加信息，只能减少」。$X$ 是原始数据，$Y$ 是中间表示，$Z$ 是它的进一步加工 —— 无论加工得多聪明，$Z$ 里关于 $X$ 的信息不会超过 $Y$。

**推论 2.2（等号条件）** 在 $X\to Y\to Z$ 下，$I(X;Y)=I(X;Z)$ 当且仅当 $X\to Z\to Y$ 也是马氏链（即 $I(X;Y\mid Z)=0$）。

**证明** 由上一条要点，$I(X;Y)-I(X;Z)=I(X;Y\mid Z)$，故等号 $\iff I(X;Y\mid Z)=0\iff X$ 与 $Y$ 在给定 $Z$ 下条件独立。$\square$

> **要点** 等号的含义是「$Z$ 已经把 $Y$ 里关于 $X$ 的信息全部装下了」，也就是统计学里的**充分性**：$Z$ 是 $X$ 关于 $Y$ 的充分统计量。数据处理不等式 + 等号刻画，就是「任何统计量都不比充分统计量更好」的精确表述。

**推论 2.3（反向也成立）** 在 $X\to Y\to Z$ 下，$I(Y;Z)\ge I(X;Z)$。

**证明** 由命题 1.2，$Z\to Y\to X$ 也是马氏链，对它用定理 2.1 得 $I(Z;X)\le I(Z;Y)$，再用对称性。$\square$

**推论 2.4（条件化减少互信息）** 在 $X\to Y\to Z$ 下，$I(X;Y)\ge I(X;Y\mid Z)$。

**证明** 由 $I(X;Y)-I(X;Z)=I(X;Y\mid Z)$ 与 $I(X;Z)\ge0$。$\square$

> **要点** 一般情形下 $I(X;Y\mid Z)\le I(X;Y)$ 是**错的**（见 [[01-熵、相对熵与互信息]] 里 $Z=X\oplus Y$ 的反例）。这里之所以成立，靠的是马氏性 —— 缺了 $X\to Y\to Z$，这条推论立刻失效。**别把条件形式当无条件结论背。**

**推论 2.5（马氏性的互信息刻画）** $I(X;Y,Z)=I(X;Y)$ 当且仅当 $X\to Y\to Z$。

**证明** $I(X;Y,Z)=I(X;Y)+I(X;Z\mid Y)$，故等式成立 $\iff I(X;Z\mid Y)=0\iff$ 给定 $Y$ 时 $X$ 与 $Z$ 条件独立。$\square$

**推论 2.6（确定性后处理）** 若 $X\to Y\to Z$ 且 $Z=f(Y)$，则 $I(X;Y)\ge I(X;f(Y))$。反之，对任意确定性函数 $f$，$X\to Y\to f(Y)$ 自动是马氏链，故 $I(X;Y)\ge I(X;f(Y))$ 无需任何条件。

**证明** 前半是定理 2.1 的特例。后半：给定 $Y$ 后 $f(Y)$ 完全确定，$I(X;f(Y)\mid Y)=0$。$\square$

**定理 2.7（四变量形）** 若 $X_1\to X_2\to X_3\to X_4$ 是马氏链，则

$$
I(X_1;X_4)\le I(X_2;X_3).
$$

**证明** 由命题 1.3（子链），$X_1\to X_2\to X_3$ 与 $X_1\to X_3\to X_4$ 都是马氏链。于是

$$
I(X_1;X_3)\le I(X_2;X_3)\quad(\text{对 }X_1\to X_2\to X_3\text{ 用定理 2.1}),
$$

$$
I(X_1;X_4)\le I(X_1;X_3)\quad(\text{对 }X_1\to X_3\to X_4\text{ 用定理 2.1}),
$$

两式串起来即得。$\square$

> **要点** 四变量形的意思是「跨越的步数越多，互信息越小」。注意它不能由一次 DPI 直接得到 —— 中间必须补一步 $I(X_1;X_4)\le I(X_1;X_3)$，而这一步用的是子链 $X_1\to X_3\to X_4$，不是原链。这正是命题 1.3 存在的理由。

## 3 Fano 不等式

**定义 3.1（错误概率）** 设 $X$ 取值于有限字母表 $\mathscr X$，$\hat X$ 是 $X$ 的一个估计（允许取值于任意集合，下面总假定 $\hat X\in\mathscr X$）。定义错误概率

$$
P_e=P(X\ne\hat X).
$$

显然 $P_e=0$ 时 $H(X\mid\hat X)=0$（$\hat X$ 完全决定 $X$）；$P_e$ 小时 $H(X\mid\hat X)$ 也该小。Fano 不等式把这句话量化。

**定理 3.2（Fano 不等式）** 设 $X,\hat X$ 取值于同一有限字母表 $\mathscr X$，则

$$
H(X\mid\hat X)\le H(P_e)+P_e\log\big(\lvert\mathscr X\rvert-1\big),
$$

其中 $H(P_e)$ 指二元熵 $H(P_e)=-P_e\log P_e-(1-P_e)\log(1-P_e)$。

**证明** 令 $E=\mathbf 1\{X\ne\hat X\}$，则 $E$ 是 Bernoulli$(P_e)$，$H(E)=H(P_e)$。由「条件不增加熵」与链式法则，

$$
H(X\mid\hat X)\le H(X,E\mid\hat X)=H(E\mid\hat X)+H(X\mid E,\hat X)\le H(E)+H(X\mid E,\hat X).
$$

第一项 $=H(P_e)$。第二项按 $E$ 的取值拆开：

$$
H(X\mid E,\hat X)=(1-P_e)\,H(X\mid E=0,\hat X)+P_e\,H(X\mid E=1,\hat X).
$$

当 $E=0$ 时 $X=\hat X$ 完全确定，该条件熵为 $0$。当 $E=1$ 时，已知 $\hat X$ 且 $X\ne\hat X$，故 $X$ 只可能在 $\mathscr X\setminus\{\hat X\}$ 中取值，这个集合只有 $\lvert\mathscr X\rvert-1$ 个元素，于是该条件熵不超过 $\log(\lvert\mathscr X\rvert-1)$。合并两项即得。$\square$

> **要点** 证明的关键是那条辅助变量 $E$：它把「$X$ 与 $\hat X$ 差多少」这个事件单独拎出来，使条件熵可以被拆成「猜错时的不确定性」与「猜对时的不确定性（为零）」两部分。**套上辅助变量再拆条件熵**，是 Fano 类不等式的标准起手式。

**推论 3.3（常用弱形式）** 在定理 3.2 的条件下，

$$
H(X\mid\hat X)\le 1+P_e\log\lvert\mathscr X\rvert .
$$

**证明** 用 $H(P_e)\le1$（二元熵的最大值）与 $\log(\lvert\mathscr X\rvert-1)\le\log\lvert\mathscr X\rvert$。$\square$

**推论 3.4（$P_e\to0$ 时条件熵趋于零）**

$$
0\le\varlimsup_{P_e\to0}H(X\mid\hat X)\le\varlimsup_{P_e\to0}H(P_e)=0 .
$$

即 $P_e\to0$ 蕴含 $H(X\mid\hat X)\to0$，而且这个收敛对字母表固定的一族 $(X,\hat X)$ 是一致的。

**证明** 对定理 3.2 的右端取上极限。字母表固定时右端只是 $P_e$ 的函数：$H(P_e)\to0$ 且 $P_e\log(\lvert\mathscr X\rvert-1)\to0$，两项都随 $P_e\to0$ 趋于 $0$。由于这个界不依赖 $(X,\hat X)$ 的具体分布，收敛对字母表固定的一族是一致。$\square$

**推论 3.5（常数猜测的界）** 设 $X\sim\{p_1,\dots,p_M\}$，$p_1=\max_i p_i$，取 $\hat X\equiv$（$p_1$ 对应的那个符号）。则

$$
H(X)\le H(p_1)+(1-p_1)\log(M-1),
$$

等号成立当且仅当 $X$ 的分布是 $\big(p_1,\frac{1-p_1}{M-1},\dots,\frac{1-p_1}{M-1}\big)$。

**证明** $\hat X\equiv x_1$ 是常数，故 $H(X\mid\hat X)=H(X)$；而 $P_e=P(X\ne x_1)=1-p_1$。此时 $E=\mathbf 1\{X\ne x_1\}$ 完全由 $X$ 决定，按 $E$ 拆分：

$$
H(X)=H(E)+H(X\mid E)=H(P_e)+(1-P_e)\cdot0+P_e\,H(X\mid E=1).
$$

$E=1$ 时 $X$ 落在 $\mathscr X\setminus\{x_1\}$ 上，只有 $M-1$ 个取值，故 $H(X\mid E=1)\le\log(M-1)$，等号成立当且仅当这 $M-1$ 个值等概率。代入即得不等式与等号条件。$\square$

**补充（逆定理要用的形状）** 推论 3.5 是一个纯熵的估计，完全不需要信息论工具就能验证：$M=2$ 时它退化为 $H(X)\le H(p_1)$，正是熵的凹性加对称性；$M$ 大时它说明「几乎确定的事件，其熵最多比 $H(p_1)$ 多一点点」。把它翻过来写成 $P_e$ 的下界，就是逆定理要卡的那个形状：

$$
P_e\ge\frac{H(X\mid\hat X)-1}{\log\lvert\mathscr X\rvert} .
$$

## 4 依概率收敛与大数定律

**定义 4.1（两种收敛）** 设 $X_n,X$ 定义在同一概率空间上。

1. **依概率收敛**：若对每个 $\varepsilon>0$ 有 $P(\lvert X_n-X\rvert\ge\varepsilon)\to0$，记 $X_n\xrightarrow{P}X$。
2. **几乎必然收敛**（依概率 $1$ 收敛）：若 $P\big(\lim_{n\to\infty}X_n=X\big)=1$，记 $X_n\xrightarrow{\text{a.s.}}X$。

> **要点** 几乎必然收敛强于依概率收敛，反之不然。典型反例是「滑动的方块」：在 $[0,1]$ 上依次取长度递减、位置滑动的小区间指示函数，它依概率收敛到 $0$，但每个点都被无穷多次覆盖，故不几乎必然收敛。本篇后面只用依概率收敛（AEP 就停在这一档）。

**定理 4.2（弱大数定律，二阶矩版）** 设 $X_1,X_2,\dots$ 独立同分布，$\mu=\mathbb E X_1$，$\sigma^2=\operatorname{Var}X_1<\infty$，记 $\bar X_n=\frac1n\sum_{i=1}^{n}X_i$。则

$$
\bar X_n\xrightarrow{P}\mu .
$$

**证明** 由独立同分布，$\mathbb E\bar X_n=\mu$，且

$$
\operatorname{Var}\bar X_n=\frac{1}{n^2}\sum_{i,j}\operatorname{Cov}(X_i,X_j)=\frac{1}{n^2}\cdot n\sigma^2=\frac{\sigma^2}{n}.
$$

对 $\bar X_n$ 用 Chebyshev 不等式：

$$
P\big(\lvert\bar X_n-\mu\rvert\ge\varepsilon\big)\le\frac{\operatorname{Var}\bar X_n}{\varepsilon^2}=\frac{\sigma^2}{n\varepsilon^2}\xrightarrow[n\to\infty]{}0 .
$$

$\square$

> **要点（只假设一阶矩的版本）** 定理 4.2 的真条件是 $\mathbb E\lvert X_1\rvert<\infty$，此时方差可能无限，Chebyshev 不够用。两种标准补法：
> 一是**截断法**：取 $Y_i=X_i\mathbf 1\{\lvert X_i\rvert\le M_n\}$，水平 $M_n$ 必须随 $n$ 一起长，且要让「某个 $\lvert X_i\rvert>M_n$」的概率 $nP(\lvert X_1\rvert>M_n)$ 与截断后的 Chebyshev 界 $M_n^2/(n\varepsilon^2)$ **同时**趋于零，再补一步 $\mathbb E Y_1=\mu_{M_n}\to\mu$（控制收敛）。配平 $M_n$ 是这条证明唯一的技术活，写错就废。
> 二是**特征函数法**：$\varphi_{\bar X_n}(t)=\varphi(t/n)^n$，由 $\mathbb E\lvert X_1\rvert<\infty$ 得 $\varphi(t)=1+i\mu t+o(t)$，于是 $\varphi_{\bar X_n}(t)\to e^{i\mu t}$，由 Lévy 连续性定理得依分布收敛到常数 $\mu$，而依分布收敛到常数就是依概率收敛。干净得多，但需要特征函数这套工具。
> **本篇两条都用不上**：AEP 里的 $i(X)=-\log p(X)$ 在有限字母表上有界，二阶矩自动有限，定理 4.2 就够。

## 5 渐近均分性（AEP）

**定理 5.1（AEP）** 设 $X_1,X_2,\dots$ 独立同分布，公共分布为 $p$，$H(X)<\infty$。则

$$
-\frac1n\log p(X_1,X_2,\dots,X_n)\xrightarrow{P}H(X).
$$

即对每个 $\varepsilon>0$，

$$
P\Big(\Big\lvert-\frac1n\log p(X_1^n)-H(X)\Big\rvert\ge\varepsilon\Big)\xrightarrow[n\to\infty]{}0 .
$$

**证明** 由独立性，

$$
-\frac1n\log p(x_1^n)=-\frac1n\sum_{i=1}^{n}\log p(x_i)=\frac1n\sum_{i=1}^{n}i(x_i),
$$

其中 $i(x)=-\log p(x)$ 是自信息。于是左端是 $n$ 个独立同分布随机变量 $i(X_1),\dots,i(X_n)$ 的算术平均。由定理 4.2，只要 $\mathbb E\lvert i(X_1)\rvert<\infty$ 就有

$$
\frac1n\sum_{i=1}^{n}i(X_i)\xrightarrow{P}\mathbb E\,i(X_1)=H(X).
$$

而有限字母表自动给出这条可积性：$\lvert\mathscr X\rvert<\infty$ 时 $i$ 有界（$\lvert i(x)\rvert=\lvert\log p(x)\rvert\le\log\lvert\mathscr X\rvert-\log p_{\min}<\infty$，其中 $p_{\min}$ 取遍正概率），有界即平方可积，定理 4.2 直接可用。$\square$

> **要点** AEP 的全部内容就是「自信息是 iid 的，所以平均自信息满足大数定律」。它没有任何深的技术，价值在于**解释**：几乎每个长序列的概率都近似等于 $2^{-nH(X)}$，概率不再是散开的，而是集中到一个数量级上。
> **补充** 有限字母表时 $H(X)<\infty$ 自动成立（熵有上界 $\log\lvert\mathscr X\rvert$），所以定理 5.1 里的这条假设只在可数字母表时才有内容。原稿把它写进条件是对的，只是不必额外验证。

## 6 $\varepsilon$-典型集

**定义 6.1（$\varepsilon$-典型集）** 设 $X_1,\dots,X_n$ 独立同分布，公共分布为 $p$。对 $\varepsilon>0$ 定义

$$
A_\varepsilon^{(n)}(X)=\Big\{x^n\in\mathscr X^n:\Big\lvert-\frac1n\log p(x^n)-H(X)\Big\rvert<\varepsilon\Big\},
$$

称其中元素为 $\varepsilon$-典型序列。

**定理 6.2（典型集的性质）** 设 $X_1,\dots,X_n$ 独立同分布，公共分布为 $p$，$H(X)<\infty$。

1. 若 $x^n\in A_\varepsilon^{(n)}(X)$，则 $2^{-n(H(X)+\varepsilon)}\le p(x^n)\le2^{-n(H(X)-\varepsilon)}$。
2. $P\big(A_\varepsilon^{(n)}(X)\big)\to1$；具体地，$n$ 充分大时 $P\big(A_\varepsilon^{(n)}(X)\big)>1-\varepsilon$。
3. $\big\lvert A_\varepsilon^{(n)}(X)\big\rvert\le2^{n(H(X)+\varepsilon)}$。
4. $n$ 充分大时 $\big\lvert A_\varepsilon^{(n)}(X)\big\rvert\ge(1-\varepsilon)2^{n(H(X)-\varepsilon)}$。

**证明** 第 1 条：典型性说的是 $\lvert-\frac1n\log p(x^n)-H\rvert<\varepsilon$，即

$$
H-\varepsilon<-\frac1n\log p(x^n)<H+\varepsilon .
$$

三边同乘 $-n$（不等号反向）得 $n(H-\varepsilon)>-\log p(x^n)>n(H+\varepsilon)$，写成 $2$ 的幂即 $2^{-n(H+\varepsilon)}<p(x^n)<2^{-n(H-\varepsilon)}$。

第 2 条：这就是定理 5.1 依概率收敛的逐字改写（把 $\ge\varepsilon$ 换成 $<\varepsilon$ 的补事件）。

第 3 条：$\sum_{x^n}p(x^n)\ge\sum_{x^n\in A}p(x^n)\ge\lvert A\rvert\cdot2^{-n(H+\varepsilon)}$，左端为 $1$，故 $\lvert A\rvert<2^{n(H+\varepsilon)}$。

第 4 条：取 $n$ 大使 $P(A)>1-\varepsilon$，则 $1-\varepsilon<P(A)\le\lvert A\rvert\cdot2^{-n(H-\varepsilon)}$，即 $\lvert A\rvert>(1-\varepsilon)2^{n(H-\varepsilon)}$。$\square$

> **要点** 1° 与 3°、4° 合起来说的是：**典型集里每个序列概率差不多，一共又只有约 $2^{nH(X)}$ 个**。所以「有效字母表」的大小是 $2^{nH(X)}$ 而不是 $\lvert\mathscr X\rvert^n$ —— 熵在这里第一次显出「数个数」的意义，而不只是「不确定度」。
> **补充（为什么 4° 需要 $n$ 充分大而 3° 不需要）** 3° 只用「概率之和不超过 1」，对每个 $n$ 都成立；4° 需要「典型集的概率接近 1」，那是 $n\to\infty$ 的结论。考试里区分这两条的时间依赖性能省掉很多无谓的常数纠缠。

## 7 信源编码定理

**定义 7.1（块编码器）** 设信源字母表为 $\mathscr X$，块长为 $n$。一个块编码器是一对映射

$$
c=(\phi,\psi),\qquad \phi:\mathscr X^n\to\mathscr M=\{0,1,\dots,M-1\},\qquad \psi:\mathscr M\to\mathscr X^n,
$$

分别称为编码函数与解码函数，码字数 $M\ge1$。定义

$$
\text{码率}\quad R(c)=\frac{\log_2 M}{n},\qquad\text{错误率}\quad P_e^{(n)}(c)=P\big(\psi(\phi(X^n))\ne X^n\big).
$$

若要求**无失真**（$P_e=0$），则 $\phi$ 必须在 $\mathscr X^n$ 上单射，此时 $M\ge\lvert\mathscr X\rvert^n$，即 $R\ge\log\lvert\mathscr X\rvert$。允许 $\varepsilon$ 的失真就能把码率压到 $H(X)$ 附近 —— 这正是下一条定理的内容。

**定义 7.2（可达码率）** 称正实数 $r$ 是可达的，若对任意 $\varepsilon>0$，存在 $n$ 与块长为 $n$ 的编码器 $c$，使得

$$
R(c)\le r+\varepsilon,\qquad P_e^{(n)}(c)\le\varepsilon .
$$

**定理 7.3（信源编码定理）**

$$
\inf\{r:r\ \text{可达}\}=H(X),
$$

并且 $H(X)$ 本身可达。

**证明** 分两部分。

**（一）可达性：$H(X)$ 可达。** 给定 $\varepsilon>0$，取 $n$ 充分大使

$$
P\big(A_\varepsilon^{(n)}\big)>1-\varepsilon,\qquad \big\lvert A_\varepsilon^{(n)}\big\rvert\le2^{n(H+\varepsilon)} .
$$

令 $M=\big\lceil2^{n(H+\varepsilon)}\big\rceil+1$，在典型集上取单射 $\phi$，把 $\mathscr X^n\setminus A_\varepsilon^{(n)}$ 的所有序列统统映到一个额外的码字上；解码 $\psi$ 在 $\phi(A_\varepsilon^{(n)})$ 上取逆，在额外码字上随便给一个序列。于是解码错误只可能在非典型序列上发生，

$$
P_e^{(n)}(c)\le P\big(\mathscr X^n\setminus A_\varepsilon^{(n)}\big)<\varepsilon,
$$

而

$$
R(c)=\frac{\log_2 M}{n}\le\frac{n(H+\varepsilon)+\log_2(1+2^{-n(H+\varepsilon)})+1}{n}\le H(X)+2\varepsilon
$$

（最后一个不等号对充分大的 $n$）。所以 $H(X)$ 可达。

**（二）逆定理。** 设有一列块长 $n\to\infty$ 的编码器 $c_n$，其错误率 $P_e^{(n)}(c_n)\to0$。则

$$
\varliminf_{n\to\infty}R(c_n)\ge H(X),
$$

特别地任何 $r<H(X)$ 都不可达。证明：固定 $n$，记 $R=R(c_n)$、$P_e=P_e^{(n)}(c_n)$，并令 $\hat X^n=\psi(\phi(X^n))$，它是 $X^n$ 的确定性函数，故 $X^n\to\hat X^n$ 是马氏链。

对 $(X^n,\hat X^n)$ 用 Fano 不等式（字母表为 $\mathscr X^n$，规模 $\lvert\mathscr X\rvert^n$）：

$$
H(X^n\mid\hat X^n)\le H(P_e)+P_e\log\big(\lvert\mathscr X\rvert^n-1\big)\le1+P_e\,n\log\lvert\mathscr X\rvert .
$$

另一方面，$\hat X^n$ 只取 $M$ 个值，故 $H(\hat X^n)\le\log M=nR$。于是

$$
nH(X)=H(X^n)=I(X^n;\hat X^n)+H(X^n\mid\hat X^n)\le H(\hat X^n)+H(X^n\mid\hat X^n)\le nR+1+P_e\,n\log\lvert\mathscr X\rvert .
$$

两边除以 $n$：

$$
H(X)\le R+\frac1n+P_e\log\lvert\mathscr X\rvert .
$$

令 $n\to\infty$。$\frac1n\to0$，$P_e\to0$ 是假设，于是右端趋于 $\varliminf R$，得

$$
H(X)\le\varliminf_{n\to\infty}R(c_n).
$$

$\square$

> **要点** 正定理（可达性）只用典型集：典型集小（$2^{n(H+\varepsilon)}$ 个），所以能编号；典型集概率大，所以错误率小。逆定理只用 Fano：Fano 说「错误率小 ⟹ $H(X^n\mid\hat X^n)$ 小」，而 $H(X^n\mid\hat X^n)$ 小意味着 $\hat X^n$ 里必须装着 $nH(X)$ 比特，装不下就必错。**两条定理各用一个工具，一个都不用另一个** —— 这个分工是信源编码这一节的结构骨架。
> **补充（与 Huffman 编码的关系）** 逐符号 Huffman 编码对任何 $p$ 都能做到 $R\le H(X)+1$，余项恒为 $1$ 比特、不随 $n$ 变小。块编码把余项换成 $O(1/n)$，代价是码表规模 $2^{nH}$ 随 $n$ 指数增长。所以实际系统里用的是「Huffman（或算术编码）逐符号 + 大块」的折中，理论上的最优只在 $n\to\infty$ 时达到。
> **补充（为什么定义里要留 $\varepsilon$ 余量）** 若把定义改成「$R(c)\le r$ 且 $P_e(c)=0$」，那么可达集是离散的（$R=\frac{\log M}{n}$ 只能取离散值），$H(X)$ 一般**不**可达 —— 例如 $X$ 均匀分布且 $H(X)=\log3$ 时，$\frac{\log M}{n}$ 永远取不到 $\log 3$。留了 $\varepsilon$ 之后下确界与最小值一致，定理才好写。

> **补充（两种「可达」定义为何等价）** 逆定理的证明用的是「码序列 + $P_e^{(n)}\to0$」这种定义（也是 Cover–Thomas 的写法），因为它能让 $n\to\infty$ 直接消掉 $\frac1n$ 项。原稿给的是「单发」定义（对每个 $\varepsilon$ 存在一个块长与编码器）。两者给出同一个门槛，理由如下。设 $r<H(X)$ 按单发定义可达，则对每个 $\varepsilon>0$ 有编码器满足 $R\le r+\varepsilon$、$P_e\le\varepsilon$，取 $\varepsilon\to0$ 得一列码。若块长无界，取 $n_k\to\infty$ 的子列，直接与上面的逆定理矛盾。若块长有界，先取子列使块长 $n$ 固定；此时 $M\le2^{n(r+1)}$ 也有界，再取子列使 $M$ 固定；块长与 $M$ 都固定时，$(\phi,\psi)$ 的组合只有有限多种，再取一次子列就能让**码本身**固定，而 $P_e\to0$ 只能逼出 $P_e=0$（有限字母表下 $P_e$ 取有限多个值）。$P_e=0$ 时 $\phi$ 在支撑上单射、$\hat X^n$ 完全恢复 $X^n$，于是 $nH(X)=H(X^n)=H(\hat X^n)\le\log M=nR$，得 $R\ge H(X)>r$，又与 $R\le r+\varepsilon\to r$ 冲突。故两种定义下不可达的集合相同。

---
## 订正与写法说明

原稿是课堂手写板书，结论全部正确，没有需要推翻的地方。下面几条是书写上要收紧的：

1. **马氏链定义里的「else 取 $0$」是多余的**：若某个 $p(x_i)=0$，则 $p(x_1^n)=0$，右边也含零因子。真正要约定的是条件概率在条件事件为零概率时的处理（不要让 $0/0$ 出现），而不是再加一个分支。
2. **反向马氏链（命题 1.2）漏了支撑条件**：严格说反向条件概率 $p(x_1\mid x_2^n)$ 只在 $p(x_2^n)>0$ 时有定义，所以「$X_n\to\cdots\to X_1$ 也是马氏链」应理解为「在 $p(x_2^n)>0$ 的点上成立」。原稿直接说「也构成 MC」，把这一步省了。
3. **Fano 不等式在 $\lvert\mathscr X\rvert=1$ 时右端出现 $\log 0$**：此时必然 $P_e=0$、$H(X\mid\hat X)=0$，按约定 $0\log0=0$ 处理即可，不构成反例，但要知道这不是「定理失效」。
4. **推论 3.4 的记号容易误读**：原稿写成「$0\le\frac{H(X\mid\hat X)}{P_e\to0}\le\frac{H(P_e)}{P_e\to0}=0$」的样子，横线上下的东西看着像分式。正确读法是**上极限**：$0\le\varlimsup_{P_e\to0}H(X\mid\hat X)\le\varlimsup_{P_e\to0}H(P_e)=0$。按分式读，量纲和结论都不对。
5. **弱大数定律的条件写的是「有限大期望」**，应为「有限期望」$\mathbb E\lvert X_1\rvert<\infty$；只假设 $\mathbb E X_1$ 存在（允许 $\pm\infty$）时结论不成立。
6. **可达码率的定义里两个 $\varepsilon$ 共用一个符号**：$R(c)\le r+\varepsilon$ 与 $P_e(c)\le\varepsilon$ 中的两个 $\varepsilon$ 其实是独立的两个任意小量。共用无碍（取小的那个即可），但证明里最好分开写，否则「$\varepsilon$ 同时控制码率与错误率」会让人以为有耦合。
7. **$\inf\{r:r\text{ 可达}\}=H(X)$ 与「$H(X)$ 可达」是两件事**，原稿并列写在一起。前者只要求下确界，后者要求等号取到 —— 后者才是正定理的完整内容，缺了它定理只能说「可以任意逼近」。
