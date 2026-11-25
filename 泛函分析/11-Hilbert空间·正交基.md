# 11 Hilbert 空间 · 正交基

> 系列：[[01-预备·Zorn引理与向量空间]] · [[02-度量空间与拓扑空间]] · [[03-拓扑向量空间]] · [[04-赋范空间]] · [[05-Banach空间]] · [[06-可分Banach空间]] · [[07-Banach空间·完备化与紧性]] · [[08-Hilbert空间·基本性质]] · [[09-Hilbert空间·最佳逼近]] · [[10-Hilbert空间·正交分解]] · [[12-连续线性算子·有界性与算子范数]] · [[13-连续线性算子·连续线性泛函与弱拓扑]] · [[14-连续线性算子·对偶空间的例子]] · [[15-连续线性算子·伴随算子]] · [[16-连续线性算子·投影算子]]
>
> 3.4 是第 3 章的收尾，也是全书第一次把 Hilbert 结构「用尽」。3.3 得到 $H=M\oplus M^\perp$，但那只拆成两块；3.4 把它拆到最细 —— 拆成一族互相垂直的**坐标轴**，且轴的数量由空间的维数决定。
>
> 三件事：
>
> 1. **正交基的存在**（定义 3.4.5、注 3.4.6）：用 Zorn 引理取极大的标准正交族。「极大」是代数条件，但它等价于分析条件 $\overline{\operatorname{span}}E=H$，这一步靠 3.3 的定理 3.3.10。
> 2. **展开定理与 Parseval 恒等式**（定理 3.4.7）：$x=\sum_\lambda\langle x,e_\lambda\rangle e_\lambda$ 且 $\lVert x\rVert^2=\sum_\lambda\lvert\langle x,e_\lambda\rangle\rvert^2$。这是勾股定理取极限。
> 3. **按维数分类**（命题 3.4.8、定理 3.4.11）：任意两组正交基等势，$\dim H$ 良定；两个 Hilbert 空间同构 $\iff$ 维数相同。于是**可分的无穷维 Hilbert 空间只有一个**，就是 $\ell^2$。

## 3.4 正交基

### 3.4.1 正交族与标准正交族

**定义 3.4.1** 设 $H$ 是 Hilbert 空间，$E=\{e_i:i\in I\}$ 是 $H$ 中的一族元素（下标集 $I$ 任意，可数或不可数）。

**(a)** $E$ 称为**正交族**（orthogonal family），若 $i\ne j$ 时 $\langle e_i,e_j\rangle=0$。

**(b)** 正交族 $E$ 称为**标准正交族**（orthonormal family），若 $\lVert e_i\rVert=1$ 对每个 $i\in I$ 成立；换言之

$$
\langle e_i,e_j\rangle=\begin{cases}1,&i=j,\\0,&i\ne j.\end{cases}
$$

由正交族造标准正交族只需逐项除以范数：$E$ 正交且每项非零时 $\{e_i/\lVert e_i\rVert\}$ 标准正交。所以后面凡涉及「正交」，都可以不失一般性地假定「标准正交」；反过来去掉单位化会省掉一堆常数，故证明里两种写法混用。

**例 3.4.2** 在 $\ell^2$ 中取

$$
e_1=(1,0,0,\dots),\qquad e_2=(0,1,0,\dots),\qquad\dots
$$

则 $\{e_n:n\in\mathbb N\}$ 是 $\ell^2$ 的标准正交族。它就是后面要证的「$\ell^2$ 的标准正交基」。

**例 3.4.3** 在实 Hilbert 空间 $L^2[-\pi,\pi]$ 中，函数族

$$
\Bigl\{\frac{1}{\sqrt{2\pi}},\ \frac{1}{\sqrt\pi}\sin t,\ \frac{1}{\sqrt\pi}\cos t,\ \frac{1}{\sqrt\pi}\sin 2t,\ \frac{1}{\sqrt\pi}\cos 2t,\ \dots\Bigr\}
$$

是标准正交族。直接算：$\int_{-\pi}^\pi\sin mt\sin nt\,dt=\pi\delta_{mn}$（$m,n\ge1$），$\int_{-\pi}^\pi\cos mt\cos nt\,dt=\pi\delta_{mn}$，$\int_{-\pi}^\pi\sin mt\cos nt\,dt=0$，$\int_{-\pi}^\pi 1\,dt=2\pi$，而 $\int_{-\pi}^\pi\cos nt\,dt=\int_{-\pi}^\pi\sin nt\,dt=0$ 保证它们都与常数 $1$ 正交。

> **要点（为什么常数项配 $\sqrt{2\pi}$ 而不是 $\sqrt\pi$）** 因为常数函数在 $[-\pi,\pi]$ 上的「长度」是 $\sqrt{2\pi}$，而 $\sin nt,\cos nt$ 的长度是 $\sqrt\pi$ —— 周期平均的意义上常数只占「半个」自由度。这个 $1/\sqrt2$ 的差别会在例 3.4.14 里变成 Fourier 级数中那个看似突兀的 $a_0/2$。

**命题 3.4.4** 若 $(e_i)_{i\in I}$ 是**非零**元素构成的正交族，则它线性无关。

**证明。** 任取有限个互不相同的下标 $i_1,\dots,i_n$ 与标量 $\lambda_1,\dots,\lambda_n$，设

$$
\lambda_1e_{i_1}+\lambda_2e_{i_2}+\cdots+\lambda_ne_{i_n}=0 .
$$

两边与 $e_{i_1}$ 作内积（内积对第一变量线性、对第二变量共轭线性，故可逐项展开）：

$$
\lambda_1\langle e_{i_1},e_{i_1}\rangle=\Bigl\langle \lambda_1e_{i_1}+\cdots+\lambda_ne_{i_n},\,e_{i_1}\Bigr\rangle=0 ,
$$

其中交叉项 $\langle e_{i_k},e_{i_1}\rangle=0$（$k\ge2$）已被消掉。由 $\langle e_{i_1},e_{i_1}\rangle=\lVert e_{i_1}\rVert^2\ne0$ 得 $\lambda_1=0$。对 $k=2,\dots,n$ 重复同一论证得 $\lambda_k=0$。故族线性无关。$\square$

> **要点（正交带来的最大好处）** 一般基下求坐标要解线性方程组；正交基下坐标就是一个内积 $c_\lambda=\langle x,e_\lambda\rangle$。命题 3.4.4 说明正交族「已经」是线性无关的，因此正交基是**最省事的那种基**：既无关，又能张成（稠密地）。

> **补充（正交族与 Hamel 基的区别）** 命题 3.4.4 只说正交族是**代数**线性无关的。一个无穷维 Hilbert 空间的正交基通常**不是** Hamel 基：$\ell^2$ 的标准基 $\{e_n\}$ 的有限线性组合只给出有限支撑序列，张不出 $e=(1,1/2,1/3,\dots)$。Hamel 基要用 Zorn 引理在全体线性无关子集上取极大，其势远大于 $\aleph_0$（对 $\ell^2$ 是 $2^{\aleph_0}$）。所以本节的「基」始终是**拓扑意义**的：$\overline{\operatorname{span}}E=H$，允许无穷级数。

### 3.4.2 正交基与存在性

**定义 3.4.5** Hilbert 空间 $H$ 的**正交基**（orthonormal basis）是指 $H$ 中一个**极大**的标准正交族 $E$：不存在标准正交族 $E'$ 使 $E\subsetneq E'$。

**注 3.4.6** 由 Zorn 引理，每个 Hilbert 空间都有正交基。事实上，把附录定理 A.4.1 的证明稍作修改即可。

> **补充（Zorn 论证的细节）** 令 $\mathcal P$ 为 $H$ 的全部标准正交族，按包含关系 $\subseteq$ 偏序。$H\ne\{0\}$ 时 $\mathcal P\ne\varnothing$（任取单位向量）。设 $\{E_i\}_{i\in I}$ 是 $\mathcal P$ 中的一条链，令 $E=\bigcup_iE_i$：任取 $e,f\in E$，存在 $i$ 使 $e,f\in E_i$（链的性质），故 $\langle e,f\rangle=\delta_{ef}$，即 $E\in\mathcal P$，且 $E$ 是链的上界。由 Zorn 引理，$\mathcal P$ 有极大元，即正交基。

极大性是一个纯代数的条件，但它等价于分析上的「张成整个空间」：

**命题（极大性的等价刻画，补充）** 设 $E$ 是 $H$ 中的标准正交族，则下列三条等价：

**(i)** $E$ 极大；

**(ii)** $\overline{\operatorname{span}}E=H$；

**(iii)** $E^\perp=\{0\}$。

**证明。** **(i)$\Rightarrow$(iii)**：设 $z\in E^\perp$，$z\ne0$。则 $E\cup\{z/\lVert z\rVert\}$ 仍是标准正交族（新元素与 $E$ 中每个元素正交，且自身单位），真包含 $E$，与极大性矛盾。故 $z=0$。

**(iii)$\Rightarrow$(ii)**：由定理 3.3.10，$E^{\perp\perp}=\overline{\operatorname{span}}E$。若 (iii) 成立，则 $\overline{\operatorname{span}}E=E^{\perp\perp}=\{0\}^\perp=H$。

**(ii)$\Rightarrow$(i)**：若 $E\subsetneq E'$ 且 $E'$ 标准正交，取 $e'\in E'\setminus E$。则 $e'\perp E$，由内积的连续性 $e'\perp\overline{\operatorname{span}}E=H$，特别 $e'\perp e'$，即 $\lVert e'\rVert^2=0$，与 $\lVert e'\rVert=1$ 矛盾。$\square$

> **要点（这条等价性是本节的枢纽）** 定义 3.4.5 用「极大」下定义，是因为它便于用 Zorn 引理证明存在性；而定理 3.4.7 的证明要用的是 (ii)。两者之间的桥梁正是 3.3 的核心定理 3.3.10（$S^{\perp\perp}=\overline{\operatorname{span}}S$）—— 这也解释了为什么 3.4 必须排在 3.3 之后。

现在可以确认两个基本例子：

- 例 3.4.2 的 $\{e_n\}$ 是 $\ell^2$ 的正交基：若 $x=(x_n)\perp e_n$ 对一切 $n$，则 $x_n=\langle x,e_n\rangle=0$，故 $x=0$，即 (iii) 成立。
- 例 3.4.3 的三角系是 $L^2[-\pi,\pi]$ 的正交基。这是分析中的标准结论（三角多项式在连续函数中一致稠密，连续函数在 $L^2$ 中稠密），例 3.4.14 会给出与例 3.4.13 平行的论证。

### 3.4.3 展开定理与 Parseval 恒等式

先单独拿出证明里反复出现的那条不等式。教材把它写成定理 3.4.7 证明中的一个 Claim，但它有独立的名字。

**引理（Bessel 不等式，补充）** 设 $E=\{e_\lambda:\lambda\in\Lambda\}$ 是 $H$ 中的标准正交族，$x\in H$。则对任意有限子集 $F\subseteq\Lambda$，

$$
\sum_{\lambda\in F}\lvert\langle x,e_\lambda\rangle\rvert^2=\Bigl\lVert\sum_{\lambda\in F}\langle x,e_\lambda\rangle e_\lambda\Bigr\rVert^2\le\lVert x\rVert^2 .
$$

并且集合 $\{\lambda\in\Lambda:\langle x,e_\lambda\rangle\ne0\}$ 至多可数。

**证明。** 记 $x_F=\sum_{\lambda\in F}\langle x,e_\lambda\rangle e_\lambda$。对 $\mu\in F$，

$$
\langle x-x_F,e_\mu\rangle=\langle x,e_\mu\rangle-\sum_{\lambda\in F}\langle x,e_\lambda\rangle\langle e_\lambda,e_\mu\rangle=\langle x,e_\mu\rangle-\langle x,e_\mu\rangle=0 ,
$$

故 $x-x_F\perp x_F$（$x_F$ 是 $\{e_\lambda\}_{\lambda\in F}$ 的线性组合）。由勾股定理（命题 3.3.3），

$$
\lVert x\rVert^2=\lVert x-x_F\rVert^2+\lVert x_F\rVert^2\ge\lVert x_F\rVert^2 ,
$$

而 $\lVert x_F\rVert^2=\sum_{\lambda\in F}\lvert\langle x,e_\lambda\rangle\rvert^2$（$F$ 内两两正交且单位）。这就得到不等式。

至于可数性：对每个 $k\in\mathbb N$，集合 $F_k=\{\lambda:\lvert\langle x,e_\lambda\rangle\rvert>1/k\}$ 满足

$$
\frac{\#F_k}{k^2}<\sum_{\lambda\in F_k}\lvert\langle x,e_\lambda\rangle\rvert^2\le\lVert x\rVert^2 ,
$$

故 $\#F_k\le k^2\lVert x\rVert^2<\infty$；而 $\{\lambda:\langle x,e_\lambda\rangle\ne0\}=\bigcup_kF_k$ 是可数多个有限集之并。$\square$

**定理 3.4.7** 设 $E=\{e_\lambda:\lambda\in\Lambda\}$ 是 Hilbert 空间 $H$ 的正交基。则每个 $x\in H$ 可以**唯一**地写成可数求和

$$
x=\sum_{\lambda\in\Lambda}c_\lambda e_\lambda ,
$$

即至多可数多个 $c_\lambda$ 非零，这里 $c_\lambda=\langle x,e_\lambda\rangle\in\mathbb K$。并且

$$
\lVert x\rVert^2=\sum_{\lambda\in\Lambda}\lvert\langle x,e_\lambda\rangle\rvert^2 .
$$

**证明。** 分四步。

**第一步（可数化）。** 由 Bessel 引理，$\Lambda_0=\{\lambda:\langle x,e_\lambda\rangle\ne0\}$ 至多可数。若 $\Lambda_0=\varnothing$，则 $x\perp E$，由命题 (iii) 得 $x=0$，定理平凡成立。故设 $\Lambda_0=\{\lambda_1,\lambda_2,\dots\}$ 是有限集或序列，记 $c_n=\langle x,e_{\lambda_n}\rangle$，$s_N=\sum_{n=1}^Nc_ne_{\lambda_n}$。

**第二步（收敛）。** 对 $N>M$，由正交性

$$
\lVert s_N-s_M\rVert^2=\Bigl\lVert\sum_{n=M+1}^Nc_ne_{\lambda_n}\Bigr\rVert^2=\sum_{n=M+1}^N\lvert c_n\rvert^2 .
$$

而 $\sum_{n=1}^\infty\lvert c_n\rvert^2\le\lVert x\rVert^2<\infty$（Bessel），故上式右端当 $M,N\to\infty$ 时趋于 $0$，即 $(s_N)$ 是 Cauchy 列。$H$ 完备，故 $s_N\to y\in H$，即 $\sum_{\lambda}c_\lambda e_\lambda$ 收敛，记其和为 $y$。

**第三步（$y=x$）。** 对任意 $\beta\in\Lambda$，由内积的连续性（定理 3.1.9）与标准正交性，

$$
\langle y,e_\beta\rangle=\lim_{N\to\infty}\langle s_N,e_\beta\rangle=\lim_{N\to\infty}\sum_{n=1}^Nc_n\langle e_{\lambda_n},e_\beta\rangle=c_\beta ,
$$

（若 $\beta\notin\Lambda_0$ 则右端为 $0$）。于是

$$
\langle x-y,e_\beta\rangle=\langle x,e_\beta\rangle-\langle y,e_\beta\rangle=0,\qquad\forall\beta\in\Lambda ,
$$

即 $x-y\in E^\perp=\{0\}$，故 $x=y$。展开式的**唯一性**也一并得到：若 $x=\sum_\lambda d_\lambda e_\lambda$，同样的计算给出 $d_\beta=\langle x,e_\beta\rangle=c_\beta$。

**第四步（Parseval）。** 由第二步的 $\lVert s_N\rVert^2=\sum_{n=1}^N\lvert c_n\rvert^2$ 与 $s_N\to x$、范数的连续性，

$$
\lVert x\rVert^2=\lim_{N\to\infty}\lVert s_N\rVert^2=\lim_{N\to\infty}\sum_{n=1}^N\lvert c_n\rvert^2=\sum_{\lambda\in\Lambda}\lvert\langle x,e_\lambda\rangle\rvert^2 . \qquad\square
$$

**补充（教材用的是「网」，本质相同）** 教材不先把 $\Lambda_0$ 可数化，而是考虑按包含关系排序的网

$$
\Bigl\{x_F=\sum_{\lambda\in F}\langle x,e_\lambda\rangle e_\lambda:\ F\subseteq\Lambda\ \text{有限}\Bigr\},
$$

用 Bessel 说明 $\lVert x_F\rVert^2\le\lVert x\rVert^2$，再用「$\lVert x\rVert^2=\lVert x-x_F\rVert^2+\lVert x_F\rVert^2$」与 $x-x_F\perp x_F$ 挤出一致收敛。两条路线等价：**先证可数、再跑 Cauchy 判据**，比处理网更容易写清楚，本文采用后者。教材的多余一步是用 $\{F_0\subset F\}$ 的网估计 $\lVert x-\sum_{\lambda\in F_0}c_\lambda e_\lambda\rVert<\eta$，其实是把「可数化 + Cauchy」重讲了一遍。

**要点（Parseval 是勾股定理的极限形式）** 有限维时它就是 $\lVert x\rVert^2=\sum_{i=1}^n\lvert x_i\rvert^2$；无穷维时它说的是「**能量守恒**」：$x$ 的总能量等于它在各坐标轴上投影的能量之和，没有泄漏。等价的写法是

$$
\bigl\lVert x-\sum_{\lambda\in F}c_\lambda e_\lambda\bigr\rVert^2=\lVert x\rVert^2-\sum_{\lambda\in F}\lvert c_\lambda\rvert^2 ,
$$

右端随 $F$ 增大单调降到 $0$，这就是「最佳逼近误差 = 总能量 − 已取能量」的一般形式（对照 [[09-Hilbert空间·最佳逼近]] 的定理 3.2.2）。

**补充（正交基的五种等价说法）** 设 $E=\{e_\lambda\}$ 标准正交，则下列等价：

1. $E$ 极大（定义 3.4.5）；
2. $\overline{\operatorname{span}}E=H$；
3. $E^\perp=\{0\}$；
4. 对一切 $x\in H$ 有 $x=\sum_\lambda\langle x,e_\lambda\rangle e_\lambda$；
5. 对一切 $x\in H$ 有 $\lVert x\rVert^2=\sum_\lambda\lvert\langle x,e_\lambda\rangle\rvert^2$。

定理 3.4.7 证的是 (1)$\Rightarrow$(4)(5)；反过来 (4)$\Rightarrow$(3)（若 $x\perp E$ 则 (4) 给 $x=0$）与 (5)$\Rightarrow$(3)（若 $x\perp E$ 则 (5) 给 $\lVert x\rVert=0$）都只有一行，再由前面的命题闭环。**关键在于「对一切 $x$」这个量词**：对单个 $x$，把 $x=x_1+x_2$（$x_1\in\overline{\operatorname{span}}E$，$x_2\perp E$）分解后，Bessel 不等式的差额正好是 $\lVert x_2\rVert^2$，故

$$
\lVert x\rVert^2=\sum_\lambda\lvert\langle x,e_\lambda\rangle\rvert^2\iff x\in\overline{\operatorname{span}}E ,
$$

所以等式可以「逐点」判断某个 $x$ 是否落在闭包里 —— 这是实际计算中检验展开是否成立最省事的办法。

### 3.4.4 维数

**命题 3.4.8** Hilbert 空间 $H$ 的任意两组正交基等势。这个公共的基数称为 $H$ 的（Hilbert 空间）**维数**，记作 $\dim H$。

**证明。** 设 $E_1,E_2$ 是两组正交基。

**有限情形。** 若 $E_1$ 有限，则由极大性 $\overline{\operatorname{span}}E_1=H$，而 $\operatorname{span}E_1$ 是有限维子空间因而闭，故 $\operatorname{span}E_1=H$，即 $\dim H=\#E_1<\infty$，且 $E_1$ 是 $H$ 的 Hamel 基。于是 $E_2$ 作为线性无关组（命题 3.4.4）满足 $\#E_2\le\dim H=\#E_1$，特别 $E_2$ 也有限；对它重复上面的论证得 $\operatorname{span}E_2=H$，故 $E_2$ 同样是 Hamel 基。两组 Hamel 基等势（有限维线性代数），于是 $\#E_1=\#E_2$。

**无限情形。** 设 $E_1,E_2$ 都无限。对每个 $e\in E_1$ 令

$$
E_e=\{f\in E_2:\langle f,e\rangle\ne0\} .
$$

把定理 3.4.7 用在正交基 $E_2$ 与向量 $e$ 上：$e=\sum_{f\in E_2}\langle e,f\rangle f$ 且 $\lVert e\rVert^2=\sum_{f\in E_2}\lvert\langle e,f\rangle\rvert^2=1<\infty$，故只有可数多个 $\langle e,f\rangle\ne0$，即 $E_e$ 至多可数。

另一方面，把定理 3.4.7 用在正交基 $E_1$ 与向量 $f\in E_2$ 上，得 $f=\sum_{e\in E_1}\langle f,e\rangle e$；因 $f\ne0$，必有某个 $e\in E_1$ 使 $\langle f,e\rangle\ne0$，即 $f\in E_e$。于是 $E_2=\bigcup_{e\in E_1}E_e$，从而

$$
\#E_2\le\sum_{e\in E_1}\#E_e\le\#E_1\cdot\aleph_0=\#E_1 ,
$$

末一步用了 $\#E_1$ 无限。由对称性同样有 $\#E_1\le\#E_2$，故 $\#E_1=\#E_2$。$\square$

> **要点（基数算术在这里只用到一条）** 无限基数 $\kappa$ 满足 $\kappa\cdot\aleph_0=\kappa$。所以「每个 $E_e$ 至多可数，而 $E_1$ 有 $\kappa$ 个」给出 $\#E_2\le\kappa$，不多不少。**这里不能反过来说「可数并至多可数」** —— 那是 $\kappa=\aleph_0$ 的特殊情形。

> **补充（教材此处写法可商榷）** 教材写的是 $\#E_2\le\sum_{e}\#E_e\le(\#E_1)^{\aleph_0}=\#E_1$。中间那步的真实理由是「$\kappa$ 个至多可数的集合之并，势 $\le\kappa\cdot\aleph_0$」，而不是「映射空间 $\kappa^{\aleph_0}$ 的大小」；两者数值上相同（$\kappa^{\aleph_0}=2^{\aleph_0\cdot\log\kappa}$ 对无限 $\kappa$ 仍等于 $\max(\kappa,2^{\aleph_0})$，一般**大于** $\kappa$），所以写成 $(\#E_1)^{\aleph_0}$ 反而会给出错误的中间量。详见「订正」第 6 条。

> **补充（$\ell^2(\Lambda)$ 与维数的例子）**
>
> - $\dim\ell^2=\aleph_0$（标准基可数）。
> - $\dim\mathbb K^n=n$；$\dim$ 有限的 Hilbert 空间必同构于 $\mathbb K^n$（定理 3.4.11）。
> - 取不可数集 $\Lambda$ 与计数测度，$H=\ell^2(\Lambda)=\{x:\Lambda\to\mathbb K:\sum_\lambda\lvert x_\lambda\rvert^2<\infty\}$，则 $\{e_\lambda\}$ 是不可数的正交基，$\dim H=\#\Lambda$。这类空间**不可分**（因为两两距离 $\sqrt2$ 的点集不可数），与 [[06-可分Banach空间]] 的判据一致。
> - 由命题 3.4.8，**可分**的 Hilbert 空间维数必 $\le\aleph_0$（见例 3.4.12(c)）。于是可分 Hilbert 空间只有 $\mathbb K^n$ 与 $\ell^2$ 两种，别无其他。

### 3.4.5 同构与酉算子

**定义 3.4.9** 两个 Hilbert 空间 $H_1,H_2$ 称为**同构**的，若存在双射线性映射 $T:H_1\to H_2$ 使

$$
\langle Tx,Ty\rangle=\langle x,y\rangle,\qquad\forall x,y\in H_1 .
$$

这样的 $T$ 称为从 $H_1$ **到** $H_2$ 的**酉算子**（unitary operator）。

**引理 3.4.10** 设 $T:H_1\to H_2$ 是 Hilbert 空间之间的**线性**映射。则 $T$ 保范数（即 $\lVert Tx\rVert=\lVert x\rVert$ 对一切 $x\in H_1$）当且仅当 $T$ 保内积（即 $\langle Tx,Ty\rangle=\langle x,y\rangle$ 对一切 $x,y$）。特别地，Hilbert 空间之间的线性双射是酉算子，当且仅当它是等距同构。

**证明。** **($\Leftarrow$)** 若保内积，取 $y=x$ 得 $\lVert Tx\rVert^2=\langle Tx,Tx\rangle=\langle x,x\rangle=\lVert x\rVert^2$。

**($\Rightarrow$)** 设保范数。由极化恒等式（复情形 (3.3)，实情形 (3.2)），对任意 $x,y$ 有

$$
4\langle Tx,Ty\rangle=\sum_{k=0}^{3}i^k\lVert Tx+i^kTy\rVert^2 .
$$

由 $T$ 线性，$Tx+i^kTy=T(x+i^ky)$；再由保范数，

$$
4\langle Tx,Ty\rangle=\sum_{k=0}^{3}i^k\lVert x+i^ky\rVert^2=4\langle x,y\rangle ,
$$

故 $\langle Tx,Ty\rangle=\langle x,y\rangle$。

最后一句：线性双射保内积 $\iff$ 保范（刚证）$\iff$ 等距同构。$\square$

> **补充（证明里两个容易被忽略的点）**
>
> 1. 极化恒等式把内积**完全**编码进范数 —— 这是「保范 $\Rightarrow$ 保内积」的唯一依据。实情形下它是 $\langle x,y\rangle=\frac14(\lVert x+y\rVert^2-\lVert x-y\rVert^2)$，复情形需要四个方向，正因为复内积有一个「相位」自由度是范数看不见的。
> 2. 第一步「$Tx+i^kTy=T(x+i^ky)$」用了 $T$ 的**线性**（不是共轭线性）。若 $T$ 只是实线性而 $H$ 是复空间，这一步就断了。

> **要点（保范只管「不发胖」，不管「不漏」）** 保范的线性映射自动是单射（$\lVert Tx\rVert=0\Rightarrow x=0$），但**不必**是满射：$\ell^2$ 上的右移 $S(x_1,x_2,\dots)=(0,x_1,x_2,\dots)$ 保范，像空间是 $\{x:x_1=0\}$。所以定义 3.4.9 要求**双射**，不能只要求保范。

**定理 3.4.11** 两个 Hilbert 空间 $H_1,H_2$ 同构，当且仅当 $\dim H_1=\dim H_2$。

**证明。** **($\Rightarrow$)** 设 $T:H_1\to H_2$ 是酉算子，$E$ 是 $H_1$ 的正交基。由引理 3.4.10，$T$ 保内积，故 $T(E)$ 仍是标准正交族；由 $T$ 是满射与线性，$\operatorname{span}T(E)=T(\operatorname{span}E)$ 在 $H_2$ 中稠密（$T$ 是等距，把稠密集映成稠密集），故 $T(E)$ 是 $H_2$ 的正交基，且 $\#T(E)=\#E$（$T$ 单射），于是 $\dim H_1=\dim H_2$。

**($\Leftarrow$)** 设 $\dim H_1=\dim H_2=\kappa$，取正交基 $E_1=\{e_\lambda:\lambda\in\Lambda\}$、$E_2=\{f_\lambda:\lambda\in\Lambda\}$（下标集相同，因为等势）。定义

$$
T\Bigl(\sum_{\lambda\in\Lambda}c_\lambda e_\lambda\Bigr)=\sum_{\lambda\in\Lambda}c_\lambda f_\lambda,\qquad \text{若}\ \sum_{\lambda\in\Lambda}c_\lambda e_\lambda\in H_1 .
$$

**良定义**：由定理 3.4.7，$x\in H_1$ 的系数 $c_\lambda=\langle x,e_\lambda\rangle$ 唯一，且 $\sum\lvert c_\lambda\rvert^2=\lVert x\rVert^2<\infty$；把定理 3.4.7 反过来用在 $H_2$ 上，$\sum_\lambda c_\lambda f_\lambda$ 收敛（同一套 Cauchy 论证，$\lVert\sum_{n=M+1}^Nc_{\lambda_n}f_{\lambda_n}\rVert^2=\sum_{n=M+1}^N\lvert c_{\lambda_n}\rvert^2$）。

**线性**：系数对 $x$ 线性且唯一。

**等距**：由两边的 Parseval，

$$
\lVert Tx\rVert^2=\sum_{\lambda}\lvert c_\lambda\rvert^2=\lVert x\rVert^2 .
$$

**满射**：任给 $y=\sum_\lambda d_\lambda f_\lambda\in H_2$（定理 3.4.7），取 $x=\sum_\lambda d_\lambda e_\lambda\in H_1$，则 $Tx=y$。

故 $T$ 是等距双射，由引理 3.4.10 是酉算子。$\square$

> **补充（教材此处一句带过）** 教材只写「由定理 3.4.7 即知 $T$ 是从 $H_1$ **到** $H_2$ 的等距」。这句话里藏了三件事：右端级数收敛（同一 Cauchy 判据）、$\lVert Tx\rVert=\lVert x\rVert$（两边 Parseval）、满射（$E_2$ 完备使每个 $y$ 都有系数）。三件事都靠定理 3.4.7，但都值得单独说一句。

**例 3.4.12**

**(a)** 每个维数 $n<\infty$ 的 Hilbert 空间 $H$ 同构于 $\mathbb K^n$。

**(b)** 每个可数无限维的 Hilbert 空间 $H$ 同构于 $\ell^2$。

**(c)** 每个可分 Hilbert 空间 $H$ 的维数可数，因而同构于 $\mathbb K^n$（$n<\infty$）或 $\ell^2$。

**(d)** $L^2[0,1]$ 同构于 $\ell^2$（因为多项式空间 $P[0,1]$ 在 $L^2[0,1]$ 中稠密，故由 (c) 得 $L^2[0,1]\cong\ell^2$）。

**证明（只证 (b)(c)，如教材）**

**(b)** 设 $(e_n)_{n\in\mathbb N}$ 是 $H$ 的正交基。对 $x\in H$，$x=\sum_{n=1}^\infty\langle x,e_n\rangle e_n$。定义

$$
\varphi:H\to\ell^2,\qquad \varphi(x)=(a_1,a_2,\dots),\quad a_n=\langle x,e_n\rangle .
$$

由 Parseval，$\lVert\varphi(x)\rVert_{\ell^2}^2=\sum_n\lvert a_n\rvert^2=\lVert x\rVert^2$，故 $\varphi$ 保范，特别 $\varphi(x)\in\ell^2$ 且 $\varphi$ 单射、连续。**满射**：任给 $\xi=(\xi_k)_{k\in\mathbb N}\in\ell^2$，令 $y_n=\sum_{k=1}^n\xi_ke_k$。对 $p\in\mathbb N$，

$$
\lVert y_{n+p}-y_n\rVert^2=\sum_{k=n+1}^{n+p}\lvert\xi_k\rvert^2\longrightarrow0\qquad(n\to\infty),
$$

因为 $\sum_k\lvert\xi_k\rvert^2<\infty$。故 $(y_n)$ 是 Cauchy 列，$H$ 完备故收敛，记极限为 $y=\sum_{k=1}^\infty\xi_ke_k$。由内积连续性 $\varphi(y)=\xi$。所以 $\varphi$ 是等距双射，$H\cong\ell^2$。

**(c)** 设 $(x_n)_{n\in\mathbb N}$ 是 $H$ 的可数稠密子集。若 $H$ 有一组不可数的正交基 $\{e_\lambda:\lambda\in\Lambda\}$，考虑开球族

$$
\bigl\{B(e_\lambda;1/2):\lambda\in\Lambda\bigr\} .
$$

任意两个不同元素的距离是 $\lVert e_\lambda-e_{\lambda'}\rVert=\sqrt2>1$，故这些球两两不交。由稠密性，每个球 $B(e_\lambda;1/2)$ 至少含一个 $x_{n_\lambda}$；由不交性，$\lambda\ne\lambda'$ 时 $n_\lambda\ne n_{\lambda'}$，即得到一个从 $\Lambda$ 到 $\mathbb N$ 的单射，与 $\Lambda$ 不可数矛盾。故维数可数；再由 (b) 与定理 3.4.11 得结论。$\square$

> **补充（(d) 里被省略的两个细节）**
>
> 1. 「$P[0,1]$ 稠密」推出可分性，但可分的定义要求**可数**稠密集。多项式全体并不可数，必须取**有理系数**多项式 $\sum_{k=0}^nq_kt^k$（$q_k\in\mathbb Q[i]$）：它可数，且由一致逼近 + $\lVert\cdot\rVert_2\le\lVert\cdot\rVert_\infty$ 仍稠密。教材这一步跳过了。
> 2. 还要 $L^2[0,1]$ 无限维。$1,t,t^2,\dots$ 线性无关（多项式次数），故维数 $=\aleph_0$。

**要点（分类定理的分量）** 定理 3.4.11 + 例 3.4.12 合起来是一个**完全分类**：可分的 Hilbert 空间在同构意义下只有

$$
\mathbb K^1,\ \mathbb K^2,\ \dots,\ \mathbb K^n,\ \dots,\ \ell^2
$$

可数多个成员。这与 Banach 空间形成尖锐对比 —— 可分的 Banach 空间有无穷多个互不同构的成员（$L^p$、$C[0,1]$、$c_0$、$\ell^1$ …），这正是 [[06-可分Banach空间]] 里「Schauder 基不唯一决定同构类」的根源：**Hilbert 空间的几何是刚性的，Banach 空间的几何不是**。

### 3.4.6 两个具体的正交基

**例 3.4.13（Legendre 多项式）** 族 $P=\{p_n(t):n=0,1,2,\dots\}$ 构成 $L^2[-1,1]$ 的正交基。$p_n$ 的构造是 Gram–Schmidt 过程的范例。

取 $A=\{1,t,t^2,\dots,t^n,\dots\}$，它在 $L^2[-1,1]$ 中线性无关但不标准正交。实空间 $L^2[-1,1]$ 的内积是

$$
\langle f,g\rangle=\int_{-1}^1f(t)g(t)\,dt .
$$

按 Gram–Schmidt 过程（练习 1）逐个正交化：$u_0=1$，$\lVert u_0\rVert=\sqrt2$，故 $p_0=1/\sqrt2$。接着 $u_1=t-\langle t,p_0\rangle p_0=t$（$\langle t,1\rangle=0$ 因 $t$ 是奇函数），$\lVert u_1\rVert=\sqrt{2/3}$，故 $p_1=\sqrt{3/2}\,t$。再往下 $u_2=t^2-\langle t^2,p_0\rangle p_0-\langle t^2,p_1\rangle p_1$，其中

$$
\langle t^2,p_0\rangle=\int_{-1}^1\frac{t^2}{\sqrt2}dt=\frac{2}{3\sqrt2},\qquad \langle t^2,p_1\rangle=\sqrt{\frac32}\int_{-1}^1t^3dt=0 ,
$$

于是 $u_2=t^2-\frac13$，$\lVert u_2\rVert^2=\int_{-1}^1\bigl(t^2-\frac13\bigr)^2dt=\frac{8}{45}$，故

$$
p_2(t)=\sqrt{\frac{45}{8}}\Bigl(t^2-\frac13\Bigr)=\sqrt{\frac52}\cdot\frac{3t^2-1}{2} .
$$

一般地

$$
p_n(t)=\sqrt{\frac{2n+1}{2}}\cdot\frac{1}{2^nn!}\frac{d^n}{dt^n}(t^2-1)^n ,
$$

括号里是 Rodrigues 公式给出的 $n$ 次 Legendre 多项式 $P_n$，故 $\lVert P_n\rVert^2=\frac{2}{2n+1}$。

**极大性。** 设 $f\in L^2[-1,1]$ 满足 $f\perp p_n$ 对一切 $n$。因 $\operatorname{span}\{1,t,t^2,\dots\}=\operatorname{span}\{p_0,p_1,p_2,\dots\}$，故

$$
\int_{-1}^1t^nf(t)\,dt=0,\qquad n=0,1,2,\dots
$$

由 Stone–Weierstrass 定理，多项式族 $P[-1,1]$ 在 $C[-1,1]$ 中一致稠密；又对 $f\in C[-1,1]$ 有 $\lVert f\rVert_2\le\sqrt2\lVert f\rVert_\infty$，故 $P[-1,1]$ 在 $C[-1,1]$ 中 $L^2$-范数稠密，而 $C[-1,1]$ 在 $L^2[-1,1]$ 中 $L^2$-范数稠密。于是 $P[-1,1]$ 在 $L^2[-1,1]$ 中稠密，从而

$$
\int_{-1}^1g(t)f(t)\,dt=0,\qquad\forall g\in L^2[-1,1] .
$$

取 $g=\bar f$（实情形即 $g=f$）得 $\int_{-1}^1\lvert f\rvert^2=0$，故 $f=0$ a.e.。这说明 $\{p_0,p_1,\dots\}$ 极大，即正交基。于是对一切 $f\in L^2[-1,1]$，

$$
f=\sum_{n=0}^\infty\langle f,p_n\rangle p_n,\qquad \lVert f\rVert^2=\sum_{n=0}^\infty\lvert\langle f,p_n\rangle\rvert^2 .
$$

**补充（Rodrigues 公式为什么正交）** 设 $m<n$，记 $Q_n(t)=(t^2-1)^n$。分部积分 $m$ 次：

$$
\int_{-1}^1t^m\frac{d^nQ_n}{dt^n}\,dt=\Bigl[t^m\frac{d^{n-1}Q_n}{dt^{n-1}}\Bigr]_{-1}^1-\int_{-1}^1mt^{m-1}\frac{d^{n-1}Q_n}{dt^{n-1}}\,dt ,
$$

边界项为零，因为 $Q_n$ 在 $\pm1$ 处有 $n$ 阶零点，故 $\frac{d^{k}Q_n}{dt^{k}}(\pm1)=0$ 对 $k<n$ 成立。重复 $m$ 次后剩下 $\pm m!\int_{-1}^1\frac{d^{n-m}Q_n}{dt^{n-m}}dt$，而 $\frac{d^{n-m}Q_n}{dt^{n-m}}$ 仍含因子 $(t^2-1)^{m}$（$n-m<n$ 阶导数的每一项都带 $(t^2-1)$ 的正幂），代入 $\pm1$ 为零。故 $\langle P_n,t^m\rangle=0$（$m<n$），即 $P_n\perp$ 一切次数更低的多项式。

**要点（Gram–Schmidt 的代价与收益）** 收益是「坐标即内积」；代价是每加一个向量要把前面所有方向都投影掉，$n$ 个向量的计算量 $O(n^2)$ 次内积。所以实际计算 Legendre 多项式用三项递推

$$
(n+1)P_{n+1}(t)=(2n+1)tP_n(t)-nP_{n-1}(t),
$$

而不是直接跑 Gram–Schmidt。这个递推可由「$tP_n$ 是 $n+1$ 次多项式，用 $P_{n+1},P_n,P_{n-1}$ 展开」得到。

**例 3.4.14（Fourier 级数）** 类似的论证说明

$$
\Bigl\{\frac{1}{\sqrt{2\pi}},\ \frac{1}{\sqrt\pi}\sin t,\ \frac{1}{\sqrt\pi}\cos t,\ \frac{1}{\sqrt\pi}\sin2t,\ \frac{1}{\sqrt\pi}\cos2t,\ \dots\Bigr\}
$$

是实 Hilbert 空间 $L^2[-\pi,\pi]$（内积 $\langle f,g\rangle=\int_{-\pi}^\pi f(t)g(t)dt$）的正交基。因此对一切 $f\in L^2[-\pi,\pi]$，

$$
f(t)=\frac{a_0}{2}+\sum_{n=1}^\infty\bigl(a_n\cos nt+b_n\sin nt\bigr) \tag{3.12}
$$

并且

$$
\lVert f\rVert^2=\pi\Bigl(\frac{a_0^2}{2}+\sum_{n=1}^\infty\bigl(a_n^2+b_n^2\bigr)\Bigr),
$$

其中

$$
a_n=\frac{1}{\pi}\int_{-\pi}^\pi f(t)\cos nt\,dt,\qquad b_n=\frac{1}{\pi}\int_{-\pi}^\pi f(t)\sin nt\,dt .
$$

数 $a_0,a_n,b_n$（$n=1,2,\dots$）称为 $f$ 的 **Fourier 系数**，(3.12) 称为 $f$ 的 **Fourier 级数展开**。

**补充（怎么从 Parseval 得到 (3.12) 与那个 $a_0/2$）** 记 $e_0=1/\sqrt{2\pi}$，$e_n=\cos nt/\sqrt\pi$，$\tilde e_n=\sin nt/\sqrt\pi$。定理 3.4.7 给

$$
f=\sum_{\text{全体}}\langle f,e\rangle e=\frac{\langle f,1\rangle}{2\pi}\cdot1+\sum_{n\ge1}\frac{\langle f,\cos nt\rangle}{\pi}\cos nt+\sum_{n\ge1}\frac{\langle f,\sin nt\rangle}{\pi}\sin nt .
$$

第一项 $\frac{\langle f,1\rangle}{2\pi}=\frac{1}{2}\cdot\frac{1}{\pi}\int f=\frac{a_0}{2}$，正是 (3.12) 的常数项；后面两项的系数正是 $a_n,b_n$。**所以 $a_0/2$ 不是约定，而是「常数方向长度 $\sqrt{2\pi}$、其余方向长度 $\sqrt\pi$」逼出来的 $1/2$。**

同理 Parseval 为

$$
\lVert f\rVert^2=\frac{\lvert\langle f,1\rangle\rvert^2}{2\pi}+\sum_{n\ge1}\frac{\lvert\langle f,\cos nt\rangle\rvert^2+\lvert\langle f,\sin nt\rangle\rvert^2}{\pi}=\pi\Bigl(\frac{a_0^2}{2}+\sum_{n\ge1}(a_n^2+b_n^2)\Bigr).
$$

**教材把这一步写成 $\lVert f\rVert^2=\frac{a_0^2}{4}+\sum(a_n^2+b_n^2)$，量纲不对**：取 $f\equiv1$，则 $a_0=2$、其余系数为零，教材给 $1$，而 $\lVert f\rVert^2=\int_{-\pi}^\pi1\,dt=2\pi$。详见「订正」第 3 条。

> **要点（三角系与 Legendre 族是同一套机器）** 两个例子的证明结构完全一样：**Gram–Schmidt 造正交族 $\to$ 说明 $\operatorname{span}$ 与某个已知稠密族相同 $\to$ 用稠密性把 $f\perp$ 整个正交族压成 $f=0$**。区别只在最后一步用哪个稠密性定理：Legendre 用 Stone–Weierstrass（多项式一致稠密），三角系用 Fejér/Weierstrass 第二定理（三角多项式一致稠密）。所以「$L^2$ 上的正交基」不是天上掉的，而是**把一致逼近定理翻译成 $L^2$ 语言**。

## 练习

**练习 1（Gram–Schmidt 过程）** 设 $\{x_n:n\in\mathbb N\}$ 是内积空间中可数线性无关族。证明存在标准正交族 $\{e_n:n\in\mathbb N\}$ 使

$$
\operatorname{span}\{e_1,\dots,e_n\}=\operatorname{span}\{x_1,\dots,x_n\},\qquad\forall n\in\mathbb N .
$$

**解。** 归纳构造。$n=1$：令 $u_1=x_1$，$e_1=u_1/\lVert u_1\rVert$（$x_1\ne0$ 由无关性），则 $\operatorname{span}\{e_1\}=\operatorname{span}\{x_1\}$。

设 $\{e_1,\dots,e_{n-1}\}$ 已造好且 $\operatorname{span}\{e_1,\dots,e_{n-1}\}=\operatorname{span}\{x_1,\dots,x_{n-1}\}$。令

$$
u_n=x_n-\sum_{i=1}^{n-1}\langle x_n,e_i\rangle e_i,\qquad e_n=\frac{u_n}{\lVert u_n\rVert} .
$$

**正交性**：对 $j<n$，

$$
\langle u_n,e_j\rangle=\langle x_n,e_j\rangle-\sum_{i=1}^{n-1}\langle x_n,e_i\rangle\langle e_i,e_j\rangle=\langle x_n,e_j\rangle-\langle x_n,e_j\rangle=0 ,
$$

故 $e_n\perp e_j$ 且 $\lVert e_n\rVert=1$。

**非退化**：若 $u_n=0$，则 $x_n=\sum_{i<n}\langle x_n,e_i\rangle e_i\in\operatorname{span}\{e_1,\dots,e_{n-1}\}=\operatorname{span}\{x_1,\dots,x_{n-1}\}$，与 $\{x_1,\dots,x_n\}$ 线性无关矛盾。故 $\lVert u_n\rVert\ne0$，$e_n$ 有定义。

**张成相同**：由定义 $e_n\in\operatorname{span}\{x_1,\dots,x_n\}$（因为 $u_n$ 是 $x_n$ 与前面 $e_i$ 的线性组合，而前面的 $e_i$ 都在 $\operatorname{span}\{x_1,\dots,x_{n-1}\}$ 中），故 $\operatorname{span}\{e_1,\dots,e_n\}\subseteq\operatorname{span}\{x_1,\dots,x_n\}$。反之 $x_n=u_n+\sum_{i<n}\langle x_n,e_i\rangle e_i=\lVert u_n\rVert e_n+\sum_{i<n}\langle x_n,e_i\rangle e_i\in\operatorname{span}\{e_1,\dots,e_n\}$，故反向包含也成立。两边维数都是 $n$，相等。$\square$

> **补充（为什么这里不需要「前缀闭」的假设）** 关键在于 $n$ 有限时 $\operatorname{span}\{x_1,\dots,x_{n-1}\}$ 是**有限维**子空间，因而闭；所以「$x_n\notin\operatorname{span}$」与「$x_n\notin\overline{\operatorname{span}}$」是一回事，$u_n$ 一定非零。练习 2 把它推到超限归纳时，这个等价就断了。

**练习 2** 设 $F=\{x_\lambda:\lambda\in\Lambda\}$ 是内积空间中（可能不可数的）线性无关族。由良序原理可设 $\Lambda$ 是良序集（附录 §A.3）。用超限归纳证明存在标准正交族 $E=\{e_\lambda:\lambda\in\Lambda\}$ 使

$$
\operatorname{span}\{e_\lambda:\lambda\in\Lambda_0\}=\operatorname{span}\{x_\lambda:\lambda\in\Lambda_0\}
$$

对 $\Lambda$ 的每个**有限**子集 $\Lambda_0$ 成立。

**解（先说结论：这题按字面是错的）。**

**反例。** 取 $H=\ell^2$，标准基记作 $(b_0,b_1,b_2,\dots)$，$\Lambda=\{1,2,\dots\}\cup\{\infty\}$（$\infty$ 为最大元），令

$$
x_n=b_n\quad(n\ge1),\qquad x_\infty=b_0+\sum_{n\ge1}\frac{b_n}{n}.
$$

$x_\infty\in\ell^2$（$\sum n^{-2}<\infty$），且这族线性无关：若 $\sum_{n\in F}c_nx_n+c_\infty x_\infty=0$（$F$ 有限），看 $b_0$ 分量得 $c_\infty=0$，再看各 $b_n$ 分量得 $c_n=0$。

设 $E=\{e_\lambda\}_{\lambda\in\Lambda}$ 满足题设。**单元集 $\Lambda_0=\{\lambda\}$ 也是有限集**，故 $\operatorname{span}\{e_\lambda\}=\operatorname{span}\{x_\lambda\}$，即两者同向：

$$
e_\lambda=\pm\frac{x_\lambda}{\lVert x_\lambda\rVert}.
$$

于是 $e_n=\pm b_n$、$e_\infty=\pm x_\infty/\lVert x_\infty\rVert$。但 $\langle x_\infty,b_1\rangle=1\ne0$，故 $e_\infty\not\perp e_1$，与 $E$ 标准正交矛盾。**所以题设的 $E$ 不存在。**

**诊断。** 题面同时要求两件事：

- $E$ 标准正交：$\langle e_\lambda,e_\mu\rangle=0$（$\lambda\ne\mu$）；
- 单元集上的张成相等：$\operatorname{span}\{e_\lambda\}=\operatorname{span}\{x_\lambda\}$，即 $e_\lambda$ 与 $x_\lambda$ 同向。

两条合起来就是「$x_\lambda$ 两两正交」。所以原族必须是正交族，此时 $e_\lambda=x_\lambda/\lVert x_\lambda\rVert$ 是唯一解，结论平凡。题面的本意应当是**闭包**意义下的张成 —— 下面给出正确版本。

**正确的命题。** 设 $F=\{x_\lambda:\lambda\in\Lambda\}$ 是内积空间中的线性无关族，$\Lambda$ 良序。则存在**子族** $\Lambda'\subseteq\Lambda$ 与标准正交族 $E=\{e_\lambda:\lambda\in\Lambda'\}$ 使

$$
\overline{\operatorname{span}}E=\overline{\operatorname{span}}F .
$$

若再把「线性无关」加强为

$$
x_\lambda\notin\overline{\operatorname{span}}\{x_\mu:\mu<\lambda\}\qquad(\forall\lambda\in\Lambda),
$$

则 $\Lambda'=\Lambda$，且对每个 $\lambda$ 有 $\overline{\operatorname{span}}\{e_\mu:\mu\le\lambda\}=\overline{\operatorname{span}}\{x_\mu:\mu\le\lambda\}$。

**证明。** 超限递归。设 $\{e_\mu:\mu<\lambda\}$ 已造好且

$$
\overline{\operatorname{span}}\{e_\mu:\mu<\lambda\}=\overline{\operatorname{span}}\{x_\mu:\mu<\lambda\},
$$

令

$$
u_\lambda=x_\lambda-\sum_{\mu<\lambda}\langle x_\lambda,e_\mu\rangle e_\mu .
$$

这个和有意义：由 Bessel 不等式 $\sum_{\mu<\lambda}\lvert\langle x_\lambda,e_\mu\rangle\rvert^2\le\lVert x_\lambda\rVert^2<\infty$，非零项至多可数，故按 Cauchy 判据收敛。

*若 $u_\lambda\ne0$*，令 $e_\lambda=u_\lambda/\lVert u_\lambda\rVert$。此时 $e_\lambda\perp e_\mu$（$\mu<\lambda$）、$\lVert e_\lambda\rVert=1$，且

$$
x_\lambda=\lVert u_\lambda\rVert e_\lambda+\sum_{\mu<\lambda}\langle x_\lambda,e_\mu\rangle e_\mu ,
$$

右端两项分别落在 $\operatorname{span}\{e_\lambda\}$ 与 $\overline{\operatorname{span}}\{e_\mu:\mu<\lambda\}$ 中，故 $x_\lambda\in\overline{\operatorname{span}}\{e_\mu:\mu\le\lambda\}$；反向由 $e_\lambda$ 的构造同理（$e_\lambda$ 是 $x_\lambda$ 与旧 $e_\mu$ 的极限组合）。于是以 $\lambda$ 结尾的闭包张成相等，归纳假设传到 $\lambda$。

*若 $u_\lambda=0$*，则 $x_\lambda\in\overline{\operatorname{span}}\{e_\mu:\mu<\lambda\}=\overline{\operatorname{span}}\{x_\mu:\mu<\lambda\}$。此时**丢掉** $x_\lambda$（取 $\lambda\notin\Lambda'$），

$$
\overline{\operatorname{span}}\{x_\mu:\mu\le\lambda\}=\overline{\operatorname{span}}\{x_\mu:\mu<\lambda\},
$$

归纳假设照样传给 $\lambda$。

*极限序数处*：$\overline{\operatorname{span}}\{x_\mu:\mu<\lambda\}=\overline{\bigcup_{\mu<\lambda}\overline{\operatorname{span}}\{x_\nu:\nu\le\mu\}}$，这一步逐项套归纳假设即得。

递归走完给出 $\overline{\operatorname{span}}E=\overline{\operatorname{span}}F$。在加强条件下 $u_\lambda=0$ 会给出 $x_\lambda\in\overline{\operatorname{span}}\{x_\mu:\mu<\lambda\}$，与假设矛盾，故 $\Lambda'=\Lambda$。$\square$

> **说明（$\Lambda'$ 为什么必须允许真包含于 $\Lambda$）** 极限序数处的缺口是实质的：那时前缀是无限集，$\operatorname{span}$ 一般**不闭**，「$x_\lambda$ 线性独立于前缀」只给出 $x_\lambda\notin\operatorname{span}$，给不出 $x_\lambda\notin\overline{\operatorname{span}}$，于是 $u_\lambda$ 可能为 $0$。取 $x_n=b_n$、$x_\infty=\sum_{n\ge1}b_n/n$ 就是这种情形：$u_\infty=0$，$e_\infty$ 根本造不出来。有限序数处没有这个问题（前缀的 $\operatorname{span}$ 是有限维、因而闭），这正是练习 1 与本题的分界。见「订正」第 7 条。

**练习 3** 设 $H$ 是可分 Hilbert 空间，证明 $H$ 中每个标准正交集可数。

**解。** 设 $E=\{e_\lambda:\lambda\in\Lambda\}$ 标准正交，$D$ 是可数稠密集。不同元素满足 $\lVert e_\lambda-e_{\lambda'}\rVert=\sqrt2$，故开球族 $\{B(e_\lambda;1/2)\}_{\lambda\in\Lambda}$ 两两不交（若 $z$ 同属两球，则 $\sqrt2\le\lVert e_\lambda-z\rVert+\lVert z-e_{\lambda'}\rVert<1$，矛盾）。每个球含 $D$ 中一点 $x_{\lambda}$；不交性保证 $\lambda\mapsto x_\lambda$ 是单射。于是 $\#\Lambda\le\#D=\aleph_0$。$\square$

> **要点** 这就是例 3.4.12(c) 的证明，只是把「正交基」放宽成「标准正交集」。注意半径必须 $<1/\sqrt2$ 的一半即 $\sqrt2/2$ 才能保证不交；取 $1/2$ 是安全的。

**练习 4** 设 $\{e_n:n\in\mathbb N\}$ 是 Hilbert 空间 $H$ 中的标准正交集。证明 $e_n\rightharpoonup0$（弱收敛），即对任意 $y\in H$ 有 $\langle e_n,y\rangle\to0$。

**解。** 固定 $y\in H$。由 Bessel 不等式，部分和 $\sum_{n=1}^N\lvert\langle y,e_n\rangle\rvert^2\le\lVert y\rVert^2$ 有界，故收敛级数

$$
\sum_{n=1}^\infty\lvert\langle y,e_n\rangle\rvert^2\le\lVert y\rVert^2<\infty ,
$$

其通项趋于零：$\lvert\langle y,e_n\rangle\rvert\to0$。于是 $\langle e_n,y\rangle=\overline{\langle y,e_n\rangle}\to0$。$\square$

> **要点（弱收敛不等于强收敛）** 对正交基 $\{e_n\}$ 而言 $\lVert e_n\rVert=1$，故 $e_n\not\to0$（强收敛意义下）。弱收敛只要求「在所有方向上投影趋于零」，$e_n$ 把质量「摊到越来越远的方向上」，这正是弱拓扑看不见的。这与 [[05-Banach空间]] 里「弱收敛 $\not\Rightarrow$ 强收敛」是同一现象。

**练习 5** 设 $(e_n)_{n\in\mathbb N}$ 是 $H$ 的正交基。证明：对任意标准正交族 $(f_n)_{n\in\mathbb N}$，若

$$
\sum_{n=1}^\infty\lVert e_n-f_n\rVert^2<\infty ,
$$

则 $(f_n)$ 也是正交基。

**解。** 记 $S=\sum_n\lVert e_n-f_n\rVert^2$。先做两个化简。

**(i)** $\lVert e_n-f_n\rVert^2=2-2\operatorname{Re}\langle e_n,f_n\rangle$，故 $S<\infty$ 等价于 $\sum_n\bigl(1-\operatorname{Re}\langle e_n,f_n\rangle\bigr)<\infty$。

**(ii)** 定义 $Vx=\sum_n\langle x,e_n\rangle f_n$。由 $(e_n)$ 是正交基与 Parseval，右端收敛且

$$
\lVert Vx\rVert^2=\sum_n\lvert\langle x,e_n\rangle\rvert^2=\lVert x\rVert^2 ,
$$

故 $V$ 是等距，$Ve_n=f_n$。令 $K=V-I$，则 $Ke_n=f_n-e_n$，

$$
\lVert K\rVert_{HS}^2:=\sum_n\lVert Ke_n\rVert^2=S .
$$

设 $M=\overline{\operatorname{span}}\{f_n\}$，$P$ 为到 $M^\perp$ 的正交投影。**要证 $M^\perp=\{0\}$。**

**小扰动情形 $S<1$。** 此时 $\lVert K\rVert\le\lVert K\rVert_{HS}=\sqrt S<1$，故 $V=I+K$ 可逆，逆由 Neumann 级数给出

$$
V^{-1}=\sum_{j=0}^\infty(-K)^j,\qquad \lVert V^{-1}\rVert\le\frac{1}{1-\sqrt S} .
$$

于是 $V$ 满射，$\operatorname{ran}V=H$，即 $\overline{\operatorname{span}}\{f_n\}=H$。由定理 3.4.7 的等价刻画，$(f_n)$ 是正交基 —— **小扰动情形至此证毕**。

**补充（$S<\infty$ 的一般情形）** 结论仍然成立，但上面的 Neumann 论证不够（$\sqrt S$ 未必 $<1$）。正确的路子是：$K$ 是 Hilbert–Schmidt 算子因而**紧**，于是 $V=I+K$ 是 Fredholm 算子且指标为 $0$；又 $V$ 等距故 $\ker V=\{0\}$，而 $\operatorname{coker}V=(\operatorname{ran}V)^\perp=M^\perp$，故

$$
\dim M^\perp=\operatorname{ind}V=0 .
$$

定量版本（把指标算成迹）是：$I-VV^*$ 是到 $M^\perp$ 的投影，故

$$
\dim M^\perp=\operatorname{tr}(I-VV^*)=-\operatorname{tr}\bigl(K+K^*+KK^*\bigr)=-\Bigl(2\operatorname{Re}\operatorname{tr}K+\lVert K\rVert_{HS}^2\Bigr),
$$

而 $\operatorname{tr}K=\sum_n\langle Ke_n,e_n\rangle=\sum_n\bigl(\langle f_n,e_n\rangle-1\bigr)=-\frac{S}{2}$，代入得 $\dim M^\perp=-(-S+S)=0$。

这条路线用到紧算子与迹类算子，超出本章范围。**若考试只考本章工具，记住 $S<1$ 的 Neumann 证明即可**；题面若真是 $S<\infty$，则要意识到它比本章工具难。见「订正」第 8 条。

**练习 6** 设 $f_n:\mathbb R\to\mathbb C$ 由

$$
f_n(x)=\pi^{-1/2}\frac{(x-i)^n}{(x+i)^{n+1}}
$$

给出。证明 $\{f_1,f_2,\dots\}$ 在 $L^2(\mathbb R)$ 中标准正交，即

$$
\int_{-\infty}^{\infty}\overline{f_m(t)}f_n(t)\,dt=\begin{cases}1,&m=n,\\0,&m\ne n.\end{cases}
$$

**解。** **可积性**：注意 $\lvert t-i\rvert=\lvert t+i\rvert=\sqrt{t^2+1}$，故

$$
\lvert f_n(t)\rvert^2=\frac{1}{\pi}\cdot\frac{(t^2+1)^n}{(t^2+1)^{n+1}}=\frac{1}{\pi(t^2+1)},
$$

与 $n$ 无关，且 $\int_{\mathbb R}\frac{dt}{\pi(t^2+1)}=1$。所以每个 $f_n\in L^2(\mathbb R)$ 且 $\lVert f_n\rVert=1$。

**正交性**：设 $k=n-m$。因 $\overline{f_m(t)}=\pi^{-1/2}\frac{(t+i)^m}{(t-i)^{m+1}}$，

$$
\overline{f_m(t)}f_n(t)=\frac{1}{\pi}\cdot\frac{(t+i)^m(t-i)^n}{(t-i)^{m+1}(t+i)^{n+1}}=\frac{1}{\pi}\cdot\frac{1}{t^2+1}\Bigl(\frac{t-i}{t+i}\Bigr)^{k} .
$$

作替换 $t=\tan\theta$，$\theta\in(-\pi/2,\pi/2)$，则 $dt=\sec^2\theta\,d\theta$，$\frac{1}{t^2+1}=\cos^2\theta$，且

$$
t+i=\frac{\sin\theta+i\cos\theta}{\cos\theta}=\frac{ie^{-i\theta}}{\cos\theta},\qquad t-i=\frac{-ie^{i\theta}}{\cos\theta},
$$

故 $\dfrac{t-i}{t+i}=-e^{2i\theta}$。代入得

$$
\int_{\mathbb R}\overline{f_m}f_n\,dt=\frac{1}{\pi}\int_{-\pi/2}^{\pi/2}(-1)^ke^{2ik\theta}\,d\theta .
$$

$k=0$ 时右端 $=\frac{1}{\pi}\cdot\pi=1$；$k\ne0$ 时

$$
\int_{-\pi/2}^{\pi/2}e^{2ik\theta}d\theta=\frac{e^{ik\pi}-e^{-ik\pi}}{2ik}=0 .
$$

故族标准正交。$\square$

> **补充（这族函数是什么）** 它们来自 Cayley 变换 $z\mapsto\frac{z-i}{z+i}$（把上半平面映到单位圆盘）的搬运：$f_n$ 是 $z^n$ 在 $L^2(\mathbb R)$ 中的对应物，构成上半平面 **Hardy 空间** $H^2(\mathbb C^+)$ 的标准正交基。注意它在 $L^2(\mathbb R)$ 中**不完整** —— 其闭包是 $H^2(\mathbb C^+)$，正交补是下半平面的 Hardy 空间，不是 $\{0\}$。所以这题只要求「标准正交」，若问「是否正交基」答案是否定的。

**练习 7** 设 $(e_n)_{n\in\mathbb N}$ 是 $H$ 中的标准正交集。证明对任意 $x\in H$ 与 $M>0$，

$$
\operatorname{card}\bigl(\{n:\lvert\langle x,e_n\rangle\rvert>M\}\bigr)\le\frac{\lVert x\rVert^2}{M^2} .
$$

**解。** 设 $F=\{n:\lvert\langle x,e_n\rangle\rvert>M\}$。对任意有限 $F_0\subseteq F$，由 Bessel 不等式

$$
\#F_0\cdot M^2<\sum_{n\in F_0}\lvert\langle x,e_n\rangle\rvert^2\le\lVert x\rVert^2 ,
$$

故 $\#F_0<\lVert x\rVert^2/M^2$ 对一切有限子集成立，取上确界得 $\#F\le\lVert x\rVert^2/M^2$。$\square$

> **要点** 这就是 Bessel 引理里「$\{\lambda:\lvert\langle x,e_\lambda\rangle\rvert>1/k\}$ 有限」的定量版，也说明**系数衰减的速度受范数控制**：$x$ 越「小」，大系数越少。Parseval 是它的极限形式（令 $M\to0$ 并把区间分成可数多层，各层计数相加得 $\sum_n\lvert\langle x,e_n\rangle\rvert^2\le\lVert x\rVert^2$）。

**练习 8** 设 $H$ 是无穷维可分 Hilbert 空间。证明 $H$ 上存在一个与原范数**不等价**的范数。

**解。** 取 $H$ 的一组正交基 $(e_n)$（可数，因可分）。**第一步：造不连续线性泛函。** 取 $H$ 的 Hamel 基 $\mathcal B$（Zorn 引理）。因 $H$ 无穷维，$\mathcal B$ 无限，故可从中取一个**可数无限**子集 $\{v_k\}_{k\ge1}$；记 $\{u_j\}_{j\in J}=\mathcal B\setminus\{v_k\}$（$J$ 任意指标集，可以空）。定义 $f$ 在基元上的值

$$
f(u_j)=0,\qquad f(v_k)=k\lVert v_k\rVert ,\quad k=1,2,\dots
$$

因 $\mathcal B$ 是 Hamel 基，$f$ 唯一地线性延拓为 $H$ 上的线性泛函。$f$ **不连续**：若连续，则存在 $C$ 使 $\lvert f(x)\rvert\le C\lVert x\rVert$ 对一切 $x$；但取 $x=v_k/\lVert v_k\rVert$ 得 $\lvert f(x)\rvert=k$，对 $k$ 无界，矛盾。

**第二步：定义新范数。** 令

$$
\lVert x\rVert'=\lVert x\rVert+\lvert f(x)\rvert .
$$

这是范数：正定性来自 $\lVert x\rVert$，齐次性与三角不等式来自 $f$ 线性；且 $\lVert x\rVert\le\lVert x\rVert'$，故 $\lVert\cdot\rVert'$ 比原范数强。

**第三步：不等价。** 若两范数等价，则存在 $C$ 使 $\lVert x\rVert'\le C\lVert x\rVert$，特别 $\lvert f(x)\rvert\le C\lVert x\rVert$，即 $f$ 连续，与第一步矛盾。故不等价。$\square$

> **补充（顺带得到一个不完备的范数）** 更强的结论：$(H,\lVert\cdot\rVert')$ 必然**不完备**。否则恒等映射 $\mathrm{id}:(H,\lVert\cdot\rVert')\to(H,\lVert\cdot\rVert)$ 是有界线性双射（$\lVert\mathrm{id}\rVert\le1$），由开映射定理（[[05-Banach空间]]）其逆有界，即两范数等价，矛盾。所以「无穷维空间上可以装一个不完备的范数」这件事，本质是**开映射定理的逆否**。

> **补充（一个更具体的构造）** 也可以用无界自伴算子：取 $(e_n)$ 与 $Tx=\sum_n n\langle x,e_n\rangle e_n$，其定义域 $D(T)=\{x:\sum_nn^2\lvert\langle x,e_n\rangle\rvert^2<\infty\}$ 在 $H$ 中稠密但不闭，$\lVert x\rVert_T=\lVert x\rVert+\lVert Tx\rVert$ 是 $D(T)$ 上的范数。但这是**子空间**上的范数，不是全空间上的；要得到全空间上的不等价范数，仍须走第一步的不连续泛函。

## 订正

1. **例 3.4.13 的「Lebesgue polynomials」应为「Legendre polynomials」**（勒让德多项式）。文中公式 $p_n(t)=\sqrt{\frac{2n+1}{2}}\frac{1}{2^nn!}\frac{d^n}{dt^n}(t^2-1)^n$ 恰是 Legendre 多项式的 Rodrigues 公式乘以归一化常数 $\sqrt{(2n+1)/2}$。Lebesgue 与 Legendre 是两个人（Henri Lebesgue / Adrien-Marie Legendre），此处纯属排印。本文已改。

2. **例 3.4.13 的 $p_0,p_1,p_2$ 与 Rodrigues 公式经核算无误**（$p_2=\sqrt{5/2}\cdot\frac{3t^2-1}{2}$，$\int_{-1}^1p_2^2=1$，与 $\lVert P_2\rVert^2=\frac25$ 相符），但教材没有说明 $p_n$ 的显式范数 $\lVert P_n\rVert^2=\frac{2}{2n+1}$ 从何而来。本文在补充里用分部积分给出了 $P_n$ 与低次多项式的正交性，配合 $\lVert P_n\rVert$ 的递推即可得该式。

3. **例 3.4.14 的 Parseval 公式量纲不对。** 教材写 $\lVert f\rVert^2=\frac{a_0^2}{4}+\sum(a_n^2+b_n^2)$，但按文中给定的内积 $\langle f,g\rangle=\int_{-\pi}^\pi fg$ 与 $a_n=\frac1\pi\int f\cos nt$，正确结果是

$$
\lVert f\rVert^2=\pi\Bigl(\frac{a_0^2}{2}+\sum_{n=1}^\infty(a_n^2+b_n^2)\Bigr).
$$

检验：$f\equiv1$ 时 $a_0=2$、其余系数为 $0$，$\lVert f\rVert^2=2\pi$；教材公式给出 $1$。进一步，教材同时写「$a_n=\langle f,\cos nt\rangle$」，与其 (3.12) 里的 $a_0/2$ 也不相容 —— 若真按 $\langle f,\cos nt\rangle$ 取系数，(3.12) 的每一项系数都要除以 $\pi$。本文按经典约定（$a_n=\frac1\pi\int f\cos nt$）统一，并给出与标准正交系的换算。

4. **引理 3.4.10 的证明里哑变量打架。** 教材写

$$
4\langle Tx,Ty\rangle=\sum_{k=0}^{3}i^k\lVert Tx+i^kTy\rVert^2=\sum_{i=0}^{3}i^k\lVert T(x+i^ky)\rVert^2=\cdots
$$

第二、三个 $\sum$ 的指标写成了 $i$，但系数仍是 $i^k$ —— 同一个式子里两个哑变量名混用。应统一为 $k$（本文已改）。另外第三步 $\lVert T(x+i^ky)\rVert=\lVert x+i^ky\rVert$ 隐含用了 $T$ 的**线性**与保范两条，教材未点明。

5. **定理 3.4.7 的证明第一段跳步。** 教材先「设 $y=\sum\langle x,e_\lambda\rangle e_\lambda$ 收敛」算出 $x-y\perp E$，然后说「由 $E$ 的极大性可得 $x=y$」。这里缺的一环是：极大性只给出 $\overline{\operatorname{span}}E=H$（命题 (ii)），还需内积的连续性才能从 $x-y\perp E$ 推出 $x-y\perp H$。教材直到后文证收敛时才引 Theorem 3.3.5，顺序上也容易让读者误以为一步到位。本文把证明重排为「可数化 $\to$ 收敛 $\to$ 展开 $\to$ Parseval」四步。

6. **命题 3.4.8 无限情形的不等式中间项写错。** 教材写 $\#E_2\le\sum_e\#E_e\le(\#E_1)^{\aleph_0}=\#E_1$。中间那步的正确理由是 $\kappa\cdot\aleph_0=\kappa$（$\kappa=\#E_1$ 无限），而非映射空间 $\kappa^{\aleph_0}$ 的大小 —— 后者一般严格大于 $\kappa$（例如 $\aleph_0^{\aleph_0}=2^{\aleph_0}>\aleph_0$）。结论不受影响，但中间量写错了对象。

7. **练习 2 按字面是错的，而且错在两处。** 第一处：$\Lambda_0$ 取单元集就得 $\operatorname{span}\{e_\lambda\}=\operatorname{span}\{x_\lambda\}$，与「$E$ 标准正交」合并等于要求原族已经两两正交 —— 反例 $x_n=b_n$、$x_\infty=b_0+\sum b_n/n$ 说明这样的 $E$ 不存在。第二处：即使把结论放宽成闭包张成，极限序数处 $u_\lambda$ 也可能为 $0$（$x_\infty=\sum b_n/n$），必须允许丢掉该元。本文已把题目改成正确的闭包版本并给出超限归纳证明。

8. **练习 5 的假设与本章工具不匹配。** 若按题面 $\sum\lVert e_n-f_n\rVert^2<\infty$，结论为真但证明需要紧算子（Fredholm 指标）或迹类算子，超出第 3 章；若把题面改成 $<1$，则 Neumann 级数即可（本文给出）。教材多半原意是后者（Paley–Wiener 型稳定性定理的标准形式）。

9. **定理 3.4.11 的「only if」方向说「plain」略快。** 「酉算子把极大标准正交族送到极大标准正交族」需要两句：保内积给出像仍是标准正交族；$T$ 满射 + 线性给出 $\operatorname{span}T(E)=T(\operatorname{span}E)$ 稠密。本文已补。

10. **例 3.4.12(d) 的「$P[0,1]$ 稠密故可数维」跳过了可分性的可数要求。** 多项式全体不可数，须取有理系数多项式才是可数稠密集；另外还需 $L^2[0,1]$ 无限维（$1,t,t^2,\dots$ 无关）。本文已补。

## 下一步

第 3 章到此结束。第 4 章从书 p73 起讲**连续线性算子**：有界性与连续性等价、算子范数、$\mathcal B(H_1,H_2)$ 的完备性、Riesz 表示定理（$H^*\cong H$）、伴随算子 $T^*$、以及自伴/正规/酉算子的谱性质。3.4 已经把「正交基」这件工具备好 —— 第 4 章会用它在 $\ell^2$ 上把算子写成（可能无限的）矩阵，并在可分情形把很多问题化归为矩阵问题。

---

系列：[[01-预备·Zorn引理与向量空间]] · [[02-度量空间与拓扑空间]] · [[03-拓扑向量空间]] · [[04-赋范空间]] · [[05-Banach空间]] · [[06-可分Banach空间]] · [[07-Banach空间·完备化与紧性]] · [[08-Hilbert空间·基本性质]] · [[09-Hilbert空间·最佳逼近]] · [[10-Hilbert空间·正交分解]] · [[12-连续线性算子·有界性与算子范数]] · [[13-连续线性算子·连续线性泛函与弱拓扑]] · [[14-连续线性算子·对偶空间的例子]] · [[15-连续线性算子·伴随算子]] · [[16-连续线性算子·投影算子]]
