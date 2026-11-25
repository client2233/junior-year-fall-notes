# 05 Banach 空间

> 系列：[[01-预备·Zorn引理与向量空间]] · [[02-度量空间与拓扑空间]] · [[03-拓扑向量空间]] · [[04-赋范空间]] · [[06-可分Banach空间]] · [[07-Banach空间·完备化与紧性]] · [[08-Hilbert空间·基本性质]] · [[09-Hilbert空间·最佳逼近]] · [[10-Hilbert空间·正交分解]] · [[11-Hilbert空间·正交基]] · [[12-连续线性算子·有界性与算子范数]] · [[13-连续线性算子·连续线性泛函与弱拓扑]] · [[14-连续线性算子·对偶空间的例子]] · [[15-连续线性算子·伴随算子]] · [[16-连续线性算子·投影算子]]
>
> 本节的核心问题只有一个：**赋范空间什么时候完备？** 答案掌握得越具体，后面（开映射定理、共鸣定理、谱理论）就越省力。

## 2.2 Banach 空间

### 2.2.1 定义

**定义 2.2.1（Cauchy 列）** 赋范空间 $E$ 中的序列 $(x_n)_{n=1}^\infty$ 称 **Cauchy 列**，若对任意 $\varepsilon>0$ 存在 $N>0$ 使

$$
\lVert x_m-x_n\rVert<\varepsilon,\qquad\text{只要}\ m,n>N.
$$

**定义 2.2.3（完备 / Banach 空间）** 赋范空间 $E$ 称**完备**的，或称为 **Banach 空间**，若 $E$ 中每个 Cauchy 列都收敛。

**例 2.2.2（完备性不是自动的）** 在赋范空间中，并非每个 Cauchy 列都收敛。例如在 $C[0,1]$ 中可找到在 $\lVert\cdot\rVert_1$ 下 Cauchy、但在 $C[0,1]$ 中不收敛的序列。取（$n\ge2$）

$$
f_n(t)=\begin{cases}
0, & 0\le t\le\dfrac12-\dfrac1n,\\[6pt]
1-n\Bigl(\dfrac12-t\Bigr), & \dfrac12-\dfrac1n\le t\le\dfrac12,\\[6pt]
1, & \dfrac12<t\le1.
\end{cases}
$$

**证明**　对 $n>m\ge2$，两个函数只可能在 $\bigl(\tfrac12-\tfrac1m,\ \tfrac12\bigr)$ 上不同（那里二者取值都在 $[0,1]$），故

$$
\lVert f_n-f_m\rVert_1=\int_0^1|f_n(t)-f_m(t)|\,dt=\int_{\frac12-\frac1m}^{\frac12}|f_n(t)-f_m(t)|\,dt\le\frac1m.
$$

所以它是 $\lVert\cdot\rVert_1$ 下的 Cauchy 列。

另一方面，设它在 $\lVert\cdot\rVert_1$ 下收敛到某连续函数 $f$。则对一切 $n>m\ge2$，

$$
\begin{aligned}
\lVert f_n-f\rVert_1&=\int_0^{\frac12}|f_n(t)-f(t)|\,dt+\int_{\frac12}^1|f(t)-1|\,dt\\
&\ge\int_0^{\frac12-\frac1n}|f(t)|\,dt+\int_{\frac12}^1|f(t)-1|\,dt.
\end{aligned}
$$

令 $n\to\infty$（左端 $\to0$），得

$$
\int_0^{\frac12-\frac1m}|f(t)|\,dt+\int_{\frac12}^1|f(t)-1|\,dt=0.
$$

故 $f(t)=0$ 当 $0\le t\le\dfrac12-\dfrac1m$，$f(t)=1$ 当 $\dfrac12\le t\le1$。这些条件对一切 $m\in\mathbb N$ 成立，故极限函数必为

$$
f(t)=\begin{cases}0,&0\le t<\dfrac12,\\[4pt] 1,&\dfrac12\le t\le1,\end{cases}
$$

而它不连续，矛盾。因此 Cauchy 列 $(f_n)_{n=2}^\infty$ 在赋范空间 $(C[0,1],\lVert\cdot\rVert_1)$ 中不收敛。$\square$

> **要点** 病根在于 $\lVert\cdot\rVert_1$ **太弱**：它看不见「单点跳变」。换成 $\lVert\cdot\rVert_\infty$ 就好了（例 2.2.4）。同一向量空间配不同范数，完备性可以完全不同 —— 这正是为什么「Banach 空间」必须写成二元组 $(E,\lVert\cdot\rVert)$。

**例 2.2.4** 例 2.2.2 中的赋范空间 $(C[0,1],\lVert\cdot\rVert_1)$ **不是** Banach 空间；但 $C[a,b]$ 配 sup 范数 $\lVert\cdot\rVert_\infty$ 是 Banach 空间。

**证明**　设 $(x_n)_{n\in\mathbb N}$ 是 $C[a,b]$ 在 sup 范数下的 Cauchy 列。对任意 $\varepsilon>0$ 存在 $N\in\mathbb N$ 使

$$
\max_{a\le t\le b}|x_{n+p}(t)-x_n(t)|<\varepsilon,\qquad\forall n>N,\ p\in\mathbb N.
$$

由底域 $\mathbb K$（$=\mathbb R$ 或 $\mathbb C$）的完备性，逐点极限

$$
x(t)=\lim_{n\to\infty}x_n(t)
$$

对每个 $t\in[a,b]$ 存在。连续函数列的**一致极限** $x$ 在 $[a,b]$ 上连续，因此 $(x_n)_{n\in\mathbb N}$ 在 $C[a,b]$ 中收敛。$\square$

> **要点** 例 2.2.4 与 2.2.2 的对比是本节的立论：**sup 范数下 Cauchy = 一致 Cauchy，一致极限保连续，所以完备**。

### 2.2.2 维数与完备性的冲突

**命题 2.2.5** 不存在具有可数（向量空间）维数的无限维 Banach 空间。

**证明**　反设无限维 Banach 空间 $E$ 有可数的（Hamel 基）$\{x_1,x_2,\dots\}$。对每个 $n\in\mathbb N$，令 $E_n$ 为 $\{x_1,\dots,x_n\}$ 的线性张成，则 $E_n$ 是 $E$ 的闭子空间且

$$
E=\bigcup_{n=1}^\infty E_n.
$$

由 **Baire 纲定理**（附录定理 C.1.3），至少有一个 $E_{n_0}$ 不是无处稠密的，即存在 $x\in E$ 与 $r>0$ 使

$$
B(x,r)\subseteq\overline{E_{n_0}}=E_{n_0}.
$$

由于 $E_{n_0}$ 是向量子空间，

$$
E=\bigcup_{k=1}^\infty kB(x,r)\subseteq\bigcup_{k=1}^\infty kE_{n_0}=E_{n_0}.
$$

于是 $\dim E=\dim E_{n_0}=n_0<\infty$，矛盾。$\square$

> **要点** 结论很有力：**Banach 空间的 Hamel 维数要么有限，要么不可数。** 换言之，用可数个基向量「拼」不出完备性 —— 完备性需要不可数多的方向。
>
> 证明是 Baire 纲定理的第一次亮相，也是本节唯一的「非构造性」论证。

**例 2.2.6** 只有有限多个非零项的序列空间 $c_{00}$ **不能配成** Banach 空间；区间 $[a,b]$ 上的多项式函数空间 $\mathcal P[a,b]$ 也不能配成 Banach 空间。

**证明**　显然 $\{e_n:n=1,2,\dots\}$ 是 $c_{00}$ 的 Hamel 基，其中 $e_n$ 是 $c_{00}$ 中第 $n$ 个坐标为 $1$、其余为 $0$ 的无穷坐标向量。故 $c_{00}$ 有可数向量空间维数，由命题 2.2.5 即得结论。多项式情形类似：$\mathcal P[a,b]$ 有可数 Hamel 基 $\{p_n(x)=x^n:n=0,1,2,\dots\}$，作为向量空间与 $c_{00}$ 同构，故同样不行。$\square$

> **补充** 注意「不能配成 Banach 空间」与「不闭」是两件事。$c_{00}\subseteq\ell^\infty$ 在既有范数下不闭（例 2.2.12），而这里说的是**无论配什么范数**都不能完备。命题 2.2.5 只涉及 Hamel 维数，与范数的选取无关。

### 2.2.3 标准 Banach 空间

**例 2.2.7** 空间 $\ell^p$ 与 $L^p[a,b]$ 在相应 $p$-范数 $\lVert\cdot\rVert_p$（$1\le p\le\infty$）下都是 Banach 空间。

**证明 (a) $1\le p<\infty$ 时 $\ell^p$ 是 Banach 空间。**

设 $x_n=(\xi_k^{(n)})_{k\in\mathbb N}$ 是 $\ell^p$ 中的 Cauchy 列。对任意 $\varepsilon>0$ 存在 $N\in\mathbb N$ 使

$$
\lVert x_m-x_n\rVert_p=\Bigl(\sum_{k=1}^\infty\bigl|\xi_k^{(n)}-\xi_k^{(m)}\bigr|^p\Bigr)^{1/p}<\varepsilon,\qquad\forall n,m\ge N.\tag{2.2}
$$

由底域的完备性，对每个固定的 $k\in\mathbb N$，Cauchy 列 $(\xi_k^{(n)})_{n\in\mathbb N}$ 收敛到某 $\xi_k\in\mathbb K$。由 (2.2)，对任意 $k_0\in\mathbb N$，

$$
\Bigl(\sum_{k=1}^{k_0}\bigl|\xi_k^{(n)}-\xi_k^{(m)}\bigr|^p\Bigr)^{1/p}<\varepsilon,\qquad\forall n,m\ge N.
$$

令 $m\to\infty$ 得

$$
\Bigl(\sum_{k=1}^{k_0}\bigl|\xi_k^{(n)}-\xi_k\bigr|^p\Bigr)^{1/p}\le\varepsilon,\qquad\forall n\ge N,\ \forall k_0\in\mathbb N.\tag{2.3}
$$

令 $x=(\xi_k)_{k\in\mathbb N}$。由 (2.3)，

$$
\Bigl(\sum_{k=1}^{k_0}|\xi_k|^p\Bigr)^{1/p}\le\Bigl(\sum_{k=1}^{k_0}\bigl|\xi_k-\xi_k^{(N)}\bigr|^p\Bigr)^{1/p}+\Bigl(\sum_{k=1}^{k_0}\bigl|\xi_k^{(N)}\bigr|^p\Bigr)^{1/p}<\epsilon+\lVert x_N\rVert_p,
$$

对 $k_0=1,2,\dots$ 成立。故 $x\in\ell^p$ 且 $\lVert x\rVert_p\le\epsilon+\lVert x_N\rVert_p<\infty$。再由 (2.3) 得 $x_n\to x$。

**(b)** 同理可证 $\ell^\infty$ 是 Banach 空间。

**(c) $1\le p<\infty$ 时 $L^p[a,b]$ 是 Banach 空间。**

设 $(f_n)_{n\in\mathbb N}$ 是 $L^p[a,b]$ 中的 Cauchy 列。对任意 $\varepsilon>0$ 存在 $N\in\mathbb N$ 使

$$
\lVert f_n-f_m\rVert_p<\varepsilon,\qquad\forall n,m\ge N.
$$

可取子列 $(f_{n_k})_{k\in\mathbb N}$ 使

$$
\lVert f_{n_{k+1}}-f_{n_k}\rVert_p<\frac1{2^k},\qquad\forall k\in\mathbb N.\tag{2.4}
$$

于是 $\sum_{k=1}^\infty\lVert f_{n_{k+1}}-f_{n_k}\rVert_p<1$。令

$$
g_m(t)=|f_{n_1}(t)|+\sum_{k=1}^m|f_{n_{k+1}}(t)-f_{n_k}(t)|,\qquad t\in[a,b].
$$

则

$$
\begin{aligned}
\lVert g_m\rVert_p&\le\lVert f_{n_1}\rVert_p+\sum_{k=1}^m\lVert f_{n_{k+1}}-f_{n_k}\rVert_p\\
&\le\lVert f_{n_1}\rVert_p+\sum_{k=1}^m\frac1{2^k}<\lVert f_{n_1}\rVert_p+1\qquad(\forall m=1,2,\dots).
\end{aligned}
$$

注意 $g_m(t)$ 正且递增，令 $g_0(t)=\lim_{m\to\infty}g_m(t)$。由 Fatou 引理，

$$
\begin{aligned}
\int_a^b|g_0|^p\,dt&=\int_a^b\liminf_{m\to\infty}|g_m|^p\,dt\le\liminf_{m\to\infty}\int_a^b|g_m|^p\,dt\\
&=\lim_{m\to\infty}\lVert g_m\rVert_p^p\le\bigl(\lVert f_{n_1}\rVert_p+1\bigr)^p.
\end{aligned}
$$

故 $g_0\in L^p[a,b]$，且

$$
g_0(t)=|f_{n_1}(t)|+\sum_{k=1}^\infty|f_{n_{k+1}}(t)-f_{n_k}(t)|
$$

在 $[a,b]$ 上 a.e. 收敛到有限正数。因此级数

$$
f(t):=f_{n_1}(t)+\sum_{k=1}^\infty\bigl(f_{n_{k+1}}(t)-f_{n_k}(t)\bigr)
$$

在 $[a,b]$ 上 a.e. 收敛到有限（实或复）数。由 $|f|\le g_0$ 得 $f\in L^p[a,b]$。又

$$
\begin{aligned}
\lVert f-f_{n_m}\rVert_p&=\Bigl\lVert f-\Bigl(f_{n_1}+\sum_{k=1}^{m-1}(f_{n_{k+1}}-f_{n_k})\Bigr)\Bigr\rVert_p
=\Bigl\lVert\sum_{k=m}^\infty(f_{n_{k+1}}-f_{n_k})\Bigr\rVert_p\\
&\le\sum_{k=m}^\infty\lVert f_{n_{k+1}}-f_{n_k}\rVert_p\le\frac1{2^{m-1}},
\end{aligned}
$$

故 $f_{n_m}\to f$。最后因 $(f_n)_{n\in\mathbb N}$ 是 Cauchy 列，得 $f_n\to f$ 于 $L^p[a,b]$。

**(d) $L^\infty[a,b]$ 是 Banach 空间。**

设 $(f_n)_{n\in\mathbb N}$ 是 $L^\infty[a,b]$ 中的 Cauchy 列。对每个 $k\in\mathbb N$ 存在 $N_k>0$ 使

$$
\lVert f_n-f_m\rVert<\frac1k,\qquad\forall n,m\ge N_k.
$$

于是可找到零测可测集 $E_k$ 使

$$
\sup_{[a,b]\setminus E_k}|f_m(t)-f_n(t)|<\frac1k,\qquad\forall n,m\ge N_k.\tag{2.5}
$$

令 $E=\bigcup_{k=1}^\infty E_k$（零测），则

$$
|f_m(t)-f_n(t)|\to0\qquad\text{关于}\ t\in[a,b]\setminus E\ \text{一致}.
$$

由 $\mathbb K$ 的完备性，极限 $f(t)=\lim_nf_n(t)$ 对一切 $t\in[a,b]\setminus E$ 存在。在 $E$ 上令 $f=0$，则 $f_n(t)\to f(t)$ 在 $[a,b]$ 上 a.e. 成立，故 $f$ 可测。在 (2.5) 中令 $m\to\infty$，

$$
\sup_{[a,b]\setminus E_k}|f(t)-f_n(t)|\le\frac1k,\qquad\forall n\ge N_k,\tag{2.6}
$$

这说明对 $n\ge N_k$ 有 $f-f_n\in L^\infty[a,b]$，从而 $f\in L^\infty[a,b]$。由 (2.6) 得

$$
\lVert f(t)-f_n(t)\rVert_\infty\le\frac1k,\qquad\forall n\ge N_k,
$$

即 $f=\lim_nf_n$ 于 $L^\infty[a,b]$。$\square$

**例 2.2.8** $c_0$ 与 $c$ 都是 $\ell^\infty$ 关于 sup 范数 $\lVert\cdot\rVert_\infty$ 的闭子空间，因此都是 Banach 空间。

> **补充（判断 Banach 空间的实用策略）** 由上两条可总结出一条省力路线：
> **$X$ 是 Banach 空间 $\Longleftarrow$ $X$ 是某已知 Banach 空间的闭子空间。**
> 验证完备性时优先找「大的」已知 Banach 空间再证闭，比直接验证 Cauchy 列收敛容易得多。

### 2.2.4 有限维情形与闭合性

**命题 2.2.9** 每个有限维赋范空间 $X$ 都完备。

**证明**　设 $\dim X=n$，作为向量空间可认为 $X=\mathbb K^n$。由定理 2.1.19，$X$ 的范数与 $2$-范数 $\lVert\cdot\rVert_2$ 等价。Cauchy 列在 $(\mathbb K^n,\lVert\cdot\rVert_2)$ 中收敛，而范数等价的不等式 (2.1) 保证同一列在 $X$ 中也收敛。故 $X$ 完备。$\square$

**命题 2.2.10** 赋范空间 $X$ 的每个有限维子空间 $M$ 都闭。

**证明**　设 $(x_n)_{n\in\mathbb N}\subseteq M$ 且 $\lim_nx_n=x\in X$。特别地 $(x_n)_{n\in\mathbb N}$ 是 Cauchy 列。$M$ 作为有限维赋范空间完备，故 $x\in M$。$\square$

> **要点** 命题 2.2.9 与 2.2.10 说明：**有限维时拓扑问题几乎消失**。无限维才是一切反例的温床 —— 下面两个例子给出「无限维子空间不闭」的标准模型。

**例 2.2.11** 设 $\mathcal P$ 是 $C[0,1]$ 中一切多项式构成的子空间。令

$$
p_n(t)=1+\frac{t}{1!}+\frac{t^2}{2!}+\dots+\frac{t^n}{n!}\in\mathcal P.
$$

由 $p_n\to e^t\notin\mathcal P$ 知无限维子空间 $\mathcal P$ 不闭。

**例 2.2.12** $c_{00}$ 在 $\ell^\infty$ 中不闭。例如取

$$
x_0=\Bigl(1,\frac12,\frac13,\dots,\frac1n,\dots\Bigr),\qquad
x_n=\Bigl(1,\frac12,\dots,\frac1n,0,0,\dots\Bigr).
$$

注意 $x_n\in c_{00}$ 对一切 $n\in\mathbb N$，且

$$
\lVert x_n-x_0\rVert_\infty=\Bigl\lVert\Bigl(0,\dots,0,\frac1{n+1},\frac1{n+2},\dots\Bigr)\Bigr\rVert_\infty=\frac1{n+1}\to0,
$$

故 $x_0\in\overline{c_{00}}\setminus c_{00}$，即 $c_{00}$ 不闭。

> **补充** 例 2.2.12 与例 2.2.6 要一起看：
> - $c_{00}$ 在 $\ell^\infty$ 中不闭 $\Rightarrow$ $(c_{00},\lVert\cdot\rVert_\infty)$ 不完备；
> - 命题 2.2.5 说 $c_{00}$（可数 Hamel 维数）**无论换成什么范数**都不能完备。
>
> 两条路殊途同归，正好互为印证。

### 2.2.5 本节速查

| 空间 | 范数 | 是否 Banach | 依据 |
|---|---|---|---|
| $\mathbb K^n$ | 任意（都等价） | ✔ | 命题 2.2.9 |
| $C[a,b]$ | $\lVert\cdot\rVert_\infty$ | ✔ | 例 2.2.4 |
| $C[0,1]$ | $\lVert\cdot\rVert_1$ | ✘ | 例 2.2.2 |
| $\ell^p$，$1\le p\le\infty$ | $\lVert\cdot\rVert_p$ | ✔ | 例 2.2.7(a)(b) |
| $L^p[a,b]$，$1\le p\le\infty$ | $\lVert\cdot\rVert_p$ | ✔ | 例 2.2.7(c)(d) |
| $c$、$c_0$ | $\lVert\cdot\rVert_\infty$ | ✔ | 例 2.2.8（$\ell^\infty$ 的闭子空间） |
| $c_{00}$ | $\lVert\cdot\rVert_\infty$ | ✘ | 例 2.2.12 不闭；例 2.2.6 可数 Hamel 维 |
| $\mathcal P[a,b]$ | $\lVert\cdot\rVert_\infty$ | ✘ | 例 2.2.11 不闭；例 2.2.6 可数 Hamel 维 |
| $BV[a,b]$ | $\lVert f\rVert=\vert f(a)\vert+V(f)$ | ✔ | 标准结果（见例 2.1.13） |

三条判据的优先级（做题时按此顺序想）：

1. **有限维？** → 立刻完备（命题 2.2.9）。
2. **能否嵌入已知 Banach 空间并证明闭？** → 例 2.2.8 的路线。
3. **若已知空间上另配一个更弱的范数** → 大概率不完备（例 2.2.2 的路线）；反之，若 Hamel 维数可数则**一定**不完备（命题 2.2.5）。

---

## 订正

1. 例 2.2.2 中 $f_n$ 的第二段原文印作「$1-n(\frac12-t)$，当 $\frac12-\frac1n\le t\le\frac12$」。该式在 $t=\frac12-\frac1n$ 处取值 $0$、在 $t=\frac12$ 处取值 $1$，与第一、三段连续衔接，无误；但教材未提连续性，需读者自行核对端点，本文在正文中已点明。
2. 例 2.2.7(c) 的尾项估计原文用 $\sum_{k=m+1}^\infty$ 配 $f_{n_m}$，索引差一位（望远镜求和 $\sum_{k=1}^m(f_{n_{k+1}}-f_{n_k})=f_{n_{m+1}}-f_{n_1}$）。本文写作尾项 $\sum_{k=m}^\infty$、界 $\le 2^{-(m-1)}$，与结论 $f_{n_m}\to f$ 一致。

---

## 下一步

后续三节已补齐：

- [[06-可分Banach空间]]（2.3）：可分的定义与例子、Schauder 基、「有基 $\Rightarrow$ 可分」及其逆不成立。
- [[07-Banach空间·完备化与紧性]]（2.4）：级数刻画完备性、完备化的存在唯一、Riesz 引理与「无限维单位球不紧」、Banach 不动点定理。
- [[08-Hilbert空间·基本性质]]（3.1）：内积的定义、Cauchy–Schwarz、平行四边形律与极化恒等式、内积连续性、直和。

第 3 章后半段（3.2 最佳逼近 → 3.3 正交分解 → 3.4 正交基）待续。

---

系列：[[01-预备·Zorn引理与向量空间]] · [[02-度量空间与拓扑空间]] · [[03-拓扑向量空间]] · [[04-赋范空间]] · [[06-可分Banach空间]] · [[07-Banach空间·完备化与紧性]] · [[08-Hilbert空间·基本性质]] · [[09-Hilbert空间·最佳逼近]] · [[10-Hilbert空间·正交分解]] · [[11-Hilbert空间·正交基]] · [[12-连续线性算子·有界性与算子范数]] · [[13-连续线性算子·连续线性泛函与弱拓扑]] · [[14-连续线性算子·对偶空间的例子]] · [[15-连续线性算子·伴随算子]] · [[16-连续线性算子·投影算子]]
