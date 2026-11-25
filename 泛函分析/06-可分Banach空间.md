# 06 可分 Banach 空间

> 系列：[[01-预备·Zorn引理与向量空间]] · [[02-度量空间与拓扑空间]] · [[03-拓扑向量空间]] · [[04-赋范空间]] · [[05-Banach空间]] · [[07-Banach空间·完备化与紧性]] · [[08-Hilbert空间·基本性质]] · [[09-Hilbert空间·最佳逼近]] · [[10-Hilbert空间·正交分解]] · [[11-Hilbert空间·正交基]] · [[12-连续线性算子·有界性与算子范数]] · [[13-连续线性算子·连续线性泛函与弱拓扑]] · [[14-连续线性算子·对偶空间的例子]] · [[15-连续线性算子·伴随算子]] · [[16-连续线性算子·投影算子]]
>
> 本节只有一条主线：**一个 Banach 空间能不能被可数集逼近？** 「可分」是把不可数断言降格为可数断言的许可证 —— 有了可数稠密集，就能把「对每个元素都成立」换成「对可数多个元素成立」。
>
> 最值得记的是一条**不对称**结论：**有 Schauder 基 $\Rightarrow$ 可分，反之不真**（Enflo 1973）。

## 2.3 可分 Banach 空间

### 2.3.1 定义

**定义 2.3.1（可分）** 赋范空间 $E$ 称**可分**的（separable），若存在 $E$ 的（范数）稠密可数子集 $M$；否则称 $E$ **不可分**（inseparable）。

> **补充（等价说法）** $\overline M=E$ 即：对每个 $x\in E$ 和每个 $\varepsilon>0$，存在 $m\in M$ 使 $\lVert x-m\rVert<\varepsilon$。等价地，每个 $x\in E$ 都是 $M$ 中某列的极限。
>
> 三个容易混淆的点：
>
> 1. **可数** 只需 $M$ 可数（有限或可数无限），**不要求** $M$ 是线性子空间，也不要求它含 $0$；
> 2. **稠密** 是相对范数拓扑讲的：同一个集合换上另一个范数，可能就不再稠密；
> 3. 可分性是**拓扑性质**（等距同构下不变），却**不是**线性性质 —— 它与范数的选取纠缠在一起（例 2.4.6(a) 会给出同一个 $C[a,b]$ 在不同范数下一完备一不完备的对照）。

### 2.3.2 例子总表

**例 2.3.2** 可分与不可分的 Banach 空间：

| 空间 | 可分性 |
|---|---|
| $c_0$、$c$ | ✔ 可分 |
| $C[a,b]$（$\lVert\cdot\rVert_\infty$） | ✔ 可分 |
| $\ell^p$，$1\le p<\infty$ | ✔ 可分 |
| $L^p[a,b]$，$1\le p<\infty$ | ✔ 可分 |
| $\ell^\infty$ | ✘ 不可分 |
| $L^\infty[a,b]$ | ✘ 不可分 |

**证明** 下面分七段给出。可分的一方是「显式构造可数稠密集」，不可分的一方用的是同一个反证模板（见 2.3.3）。

---

**(a) $c$ 可分。**

设 $x=(\xi_n)_{n\in\mathbb N}\in c$，即 $\xi_n\to\xi_0$。对任意 $\varepsilon>0$，存在 $N\in\mathbb N$ 使

$$
|\xi_n-\xi_0|<\frac\varepsilon2,\qquad\forall n>N.
$$

取有理数 $r_0,r_1,\dots,r_N\in\mathbb Q$ 使

$$
|\xi_n-r_n|<\frac\varepsilon2,\qquad n=0,1,\dots,N.
$$

令 $y=(r_1,\dots,r_N,r_0,r_0,\dots)$，即「前 $N$ 位用各自的有理数、之后全部用极限值的有理数」。则

$$
\lVert x-y\rVert_\infty=\max\Bigl\{\max_{1\le n\le N}|\xi_n-r_n|,\ \sup_{n>N}|\xi_n-r_0|\Bigr\}<\max\Bigl\{\frac\varepsilon2,\ \sup_{n>N}\bigl(|\xi_n-\xi_0|+|\xi_0-r_0|\bigr)\Bigr\}<\varepsilon.
$$

因此可数集

$$
\{(r_1,\dots,r_N,r_0,r_0,\dots):\ r_0,r_1,\dots,r_N\in\mathbb Q,\ N\in\mathbb N\}
$$

在 $c$ 中稠密，故 $c$ 可分。

**(b) $c_0$ 可分。**

$c_0$ 是 $c$ 的子空间，而 $c$ 可分，故 $c_0$ 也可分。

> **补充（为什么子空间继承可分性）** 设 $A\subseteq E$（$E$ 为可分度量空间），$M$ 是 $E$ 的可数稠密集。令

$$
M'=\Bigl\{m\in M:\ B\bigl(m,\tfrac1n\bigr)\cap A\ne\emptyset\ \text{对某个}\ n\in\mathbb N\Bigr\},
$$

> 对每个「有交」的 $(m,n)$ 从交里取一点，所得集合 $M''$ 可数且在 $A$ 中稠密（对 $a\in A$ 与 $\varepsilon>0$，取 $n>2/\varepsilon$ 与 $m\in M$ 使 $\lVert a-m\rVert<1/(2n)$，则 $B(m,1/n)\cap A\ne\emptyset$，其中选出的点与 $a$ 的距离 $<2/n<\varepsilon$）。$\square$
>
> 注意：**不可分性不继承** —— 不可分空间的子空间完全可以可分。

**(c) $C[a,b]$ 可分。**

由 **Weierstrass 逼近定理**，有理系数多项式集

$$
\Bigl\{\sum_{k=0}^{k_0}r_kt^k:\ r_0,r_1,\dots,r_{k_0}\in\mathbb Q,\ k_0\in\mathbb N\Bigr\}
$$

在 $C[a,b]$（配 sup 范数）中稠密。它是可数集，故 $C[a,b]$ 可分。

> **补充** 「有理系数多项式」之所以够用，是因为 sup 范数下可以「两步走」：先用实（复）系数多项式逼近到 $\varepsilon/2$（Weierstrass），再把每个系数换成差得足够小的有理数，系数扰动引起的一致误差不超过系数误差之和。因此关键在于**多项式集本身在 $C[a,b]$ 中稠密**，而 Weierstrass 定理提供的正是这一条。

**(d) $1\le p<\infty$ 时 $\ell^p$ 可分。**

设 $x_0=(\xi_k)_{k\in\mathbb N}\in\ell^p$。对任意 $\varepsilon>0$，存在 $k_0\in\mathbb N$ 使

$$
\sum_{k=k_0+1}^\infty|\xi_k|^p<\frac{\varepsilon^p}{2}.
$$

另一方面，存在有理数 $r_1,\dots,r_{k_0}\in\mathbb Q$ 使

$$
\sum_{k=1}^{k_0}|\xi_k-r_k|^p<\frac{\varepsilon^p}{2}.
$$

令 $y_0=(r_1,\dots,r_{k_0},0,0,\dots)$，则

$$
\lVert x_0-y_0\rVert_p^p=\sum_{k=1}^{k_0}|\xi_k-r_k|^p+\sum_{k=k_0+1}^\infty|\xi_k|^p<\varepsilon^p.
$$

于是可数集

$$
\{(r_1,\dots,r_{k_0},0,0,\dots):\ r_1,\dots,r_{k_0}\in\mathbb Q,\ k_0\in\mathbb N\}
$$

在 $\ell^p$ 中稠密，故 $\ell^p$ 可分。

> **要点** 证明是**两段拼**：尾巴用级数收敛性掐掉（这是 $\ell^p$ 的「可数维数」本质），前段用有理数逼近。$\ell^\infty$ 缺的正是第一段 —— 它的「尾巴」不受 $\lVert\cdot\rVert_\infty$ 控制，这就是它不可分的技术根源。

**(e) $1\le p<\infty$ 时 $L^p[a,b]$ 可分。**

分四步逼近，每一步都换成「更规矩」的函数。设 $f\in L^p[a,b]$，$\varepsilon>0$。

**第 1 步（截断到有界）。** 因 $|f|^p$ 可积，存在 $\delta>0$ 使

$$
\int_E|f|^p\,\mathrm dt<\Bigl(\frac\varepsilon3\Bigr)^p,\qquad\text{只要}\ \mu(E)<\delta. \tag{2.7}
$$

又 $|f|^p$（从而 $f$）在 $[a,b]$ 上 a.e. 有限，故存在 $n_0\in\mathbb N$ 使

$$
\mu\bigl(\{t\in[a,b]:|f(t)|>n_0\}\bigr)<\delta. \tag{2.8}
$$

令

$$
f_{n_0}(t)=\begin{cases}f(t),&|f(t)|\le n_0,\\[2pt] 0,&|f(t)|>n_0.\end{cases}
$$

则 $f_{n_0}\in L^p[a,b]$ 且有界，并且由 (2.7)(2.8)

$$
\lVert f-f_{n_0}\rVert_p=\Bigl(\int_{\{t:|f(t)|>n_0\}}|f|^p\,\mathrm dt\Bigr)^{1/p}<\frac\varepsilon3. \tag{2.9}
$$

**第 2 步（换成连续函数）。** 对上述有界可测函数 $f_{n_0}$，存在连续函数 $g$ 使 $|g|\le n_0$ 且

$$
\mu\bigl(\{t\in[a,b]:g(t)\ne f_{n_0}(t)\}\bigr)<\Bigl(\frac{\varepsilon}{6n_0}\Bigr)^p.
$$

于是

$$
\lVert f_{n_0}-g\rVert_p=\Bigl(\int_{\{t:g(t)\ne f_{n_0}(t)\}}|f_{n_0}-g|^p\,\mathrm dt\Bigr)^{1/p}\le\bigl[(2n_0)^p\,\mu(\{t:g(t)\ne f_{n_0}(t)\})\bigr]^{1/p}<\frac\varepsilon3. \tag{2.10}
$$

**第 3 步（换成多项式）。** 对连续函数 $g$，由 Weierstrass 定理存在有理系数多项式 $h$ 使

$$
\lVert g-h\rVert_p<\frac\varepsilon3. \tag{2.11}
$$

（此处用 $\lVert\cdot\rVert_p=(b-a)^{1/p}\lVert\cdot\rVert_\infty$ 一类的范数比较，把一致逼近换成 $L^p$ 逼近。）

**第 4 步（合并）。** 由 (2.9)–(2.11) 与三角不等式，

$$
\lVert f-h\rVert_p\le\lVert f-f_{n_0}\rVert_p+\lVert f_{n_0}-g\rVert_p+\lVert g-h\rVert_p<\varepsilon.
$$

有理系数多项式集可数，故 $L^p[a,b]$ 可分。

> **要点** 四步链是「可测 $\to$ 有界 $\to$ 连续 $\to$ 多项式」的逐级正则化，分给每级 $\varepsilon/3$。这条链在实分析里反复出现（Lusin 定理、$C_c$ 在 $L^p$ 中稠密），值得整条记下来。

**(f) $\ell^\infty$ 不可分。**

取

$$
A=\{( \xi_n)_{n\in\mathbb N}\in\ell^\infty:\ \xi_n=0\ \text{或}\ 1,\ n=1,2,\dots\}.
$$

若 $x,y\in A$ 且 $x\ne y$，则存在某坐标处二者一为 $0$、一为 $1$，故

$$
\lVert x-y\rVert_\infty=1.
$$

$A$ 与 $\{0,1\}^{\mathbb N}$ 一一对应，是**不可数**集。于是 $\{B(x,1/3):x\in A\}$ 是 $\ell^\infty$ 中一族**两两不交**的开球，且这族本身不可数。

若 $\ell^\infty$ 有可数稠密集 $M$，则每个 $B(x,1/3)$ 必含 $M$ 中的点（$x$ 到 $M$ 的距离 $<1/3$），而球两两不交，故 $M$ 至少与 $A$ 一样多 —— 不可数，矛盾。因此 $\ell^\infty$ 不可分。

**(g) $L^\infty[a,b]$ 不可分。**

取闭子区间的特征函数族

$$
B=\{\mathbf 1_{[a,t]}:\ t\in[a,b]\},
$$

它不可数。对任意 $t<s$，

$$
\lVert\mathbf 1_{[a,t]}-\mathbf 1_{[a,s]}\rVert_\infty=\lVert\mathbf 1_{(t,s]}\rVert_\infty=1.
$$

于是与 (f) 完全同样的论证说明 $L^\infty[a,b]$ 不可分。$\square$

### 2.3.3 不可分性的证法模板

> **补充（记这一条就够了）** (f)(g) 是同一个模板，凡证「$E$ 不可分」都可照搬：
>
> 1. 在 $E$ 中造一个**不可数**族 $A$；
> 2. 证明 $A$ 中两两距离有**正的下界** $\delta>0$；
> 3. 于是半径 $\delta/3$ 的球两两不交，且每个球都「必须」含稠密集中的点；
> 4. 稠密集的基数 $\ge$ $A$ 的基数 $=$ 不可数，矛盾。
>
> 反过来，可分性的构造几乎总是「**先掐尾巴，再用有理数/多项式逼近**」——(a)(c)(d)(e) 全是这个套路。
>
> 一个常用推论：**可分空间的任何连续像（商空间）可分**；而**不可分空间的子空间可以可分**（如 $\ell^\infty$ 的子空间 $c_0$）。

### 2.3.4 Schauder 基

**定义 2.3.3（Schauder 基）** 设 $(e_n)_{n\in\mathbb N}$ 是 Banach 空间 $E$ 中一列**单位向量**（$\lVert e_n\rVert=1$）。若对每个 $x\in E$ 都存在**唯一**的纯量组 $(\alpha_n)$ 使

$$
\lim_{n\to\infty}\Bigl\lVert x-\sum_{i=1}^n\alpha_ie_i\Bigr\rVert=0,
$$

即

$$
x=\sum_{n=1}^\infty\alpha_ne_n,
$$

则称 $(e_n)_{n\in\mathbb N}$ 是 $E$ 的 **Schauder 基**。显然 Schauder 基必线性无关。

> **补充（为什么「显」然线性无关）** 设有限组合 $\sum_{i=1}^{k}\lambda_ie_i=0$ 且 $\lambda_k\ne0$。另一方面 $0$ 也可写成 $\sum_{i=1}^\infty 0\cdot e_i$（部分和恒为 $0$，当然收敛到 $0$）。由**唯一性**两个表示相同，故全部 $\lambda_i=0$，矛盾。$\square$
>
> Schauder 基与 Hamel 基的关键差别：Hamel 基用**有限**线性组合，Schauder 基用**级数**（极限）。这正是为什么 Schauder 基要与范数拓扑配套，而 Hamel 基纯是代数对象（[[01-预备·Zorn引理与向量空间]]）。

**定理 2.3.4** 有 Schauder 基 $(e_n)_{n\in\mathbb N}$ 的 Banach 空间 $E$ 必可分。

**证明** 集合

$$
M=\Bigl\{\sum_{i=1}^n r_ie_i:\ r_1,\dots,r_n\in\mathbb Q,\ n\in\mathbb N\Bigr\}
$$

可数。设 $x\in E$、$\varepsilon>0$。由 Schauder 基的定义，存在 $n$ 与 $\alpha_1,\dots,\alpha_n\in\mathbb K$ 使

$$
\Bigl\lVert x-\sum_{i=1}^n\alpha_ie_i\Bigr\rVert<\frac\varepsilon2.
$$

对每个 $i$ 取 $r_i\in\mathbb Q$ 使 $|\alpha_i-r_i|<\dfrac{\varepsilon}{2n}$（复情形取 $r_i\in\mathbb Q+i\mathbb Q$）。于是

$$
\Bigl\lVert\sum_{i=1}^n\alpha_ie_i-\sum_{i=1}^nr_ie_i\Bigr\rVert\le\sum_{i=1}^n|\alpha_i-r_i|\,\lVert e_i\rVert<\frac\varepsilon2.
$$

合并得 $\bigl\lVert x-\sum_{i=1}^nr_ie_i\bigr\rVert<\varepsilon$，故 $M$ 在 $E$ 中稠密，$E$ 可分。$\square$

> **要点** 证明里 $\lVert e_i\rVert=1$ 那一项是关键：它把「系数误差」直接翻译成「向量误差」。若基向量长度不受控，$M$ 就可能不稠密。

### 2.3.5 基与可分性的不对称

> **补充（定理 2.3.4 的逆不成立）**
>
> - **可分的 Banach 空间未必有 Schauder 基。** 这是 Enflo 于 1973 年给出的反例（后来有更简洁的构造），它同时解决了 Banach 提出的「基问题」，是 20 世纪泛函分析的重要事件。
> - 常见的**特例对照**：$c_0$、$c$、$\ell^p\ (1\le p<\infty)$、$C[0,1]$ 都有 Schauder 基（见练习）；而 $\ell^\infty$ **没有** —— 但它不可分，所以卡在定理 2.3.4 的前提上，并不构成上述「可分却无基」的反例。
> - 一句话记忆：**Schauder 基是可分性的充分条件，不是必要条件。**
>
> 另外注意「$E$ 有 Schauder 基」与「$E$ 与某个经典序列空间同构」是两回事；有基的空间类比同构类宽得多。

### 2.3.6 练习

1. 设 $e_k$ 是第 $k$ 个坐标为 $1$、其余为 $0$ 的纯量列。证明 $\{e_k:k=1,2,\dots\}$ 是 Banach 空间 $(c_0,\lVert\cdot\rVert_\infty)$ 的 Schauder 基。
2. 设 $e_k$ 同上，$e=(1,1,\dots)$ 为全坐标 $1$ 的列。证明 $\{e_k:k=1,2,\dots\}\cup\{e\}$ 是 Banach 空间 $(c,\lVert\cdot\rVert_\infty)$ 的 Schauder 基。
3. 设 $e_k$ 同上。证明 $1\le p<\infty$ 时 $\{e_k:k=1,2,\dots\}$ 是 Banach 空间 $(\ell^p,\lVert\cdot\rVert_p)$ 的 Schauder 基。
4. Banach 空间 $(\ell^\infty,\lVert\cdot\rVert_\infty)$ 有 Schauder 基吗？
5. Banach 空间 $(C[0,1],\lVert\cdot\rVert_\infty)$ 有 Schauder 基。提示：定义 $\varphi_0(x)=1$、$\varphi_1(x)=x$，并对 $n=2^k+i$（$i=1,2,\dots,2^k$，$k=1,2,\dots$）令

   $$
   \varphi_n(x)=\varphi_{2^k+i}(x)=\max\Bigl\{1-2^{k+1}\Bigl|x-\frac{2i-1}{2^{k+1}}\Bigr|,\ 0\Bigr\}.
   $$

   注意 $\varphi_{2^k+i}$ 在区间 $\bigl(\frac{i-1}{2^k},\frac{i}{2^k}\bigr)$ 外恒为 $0$。对 $f\in C[0,1]$，令 $a_0(f)=f(0)$、$a_1(f)=f(1)-f(0)$，并令

   $$
   a_n(f)=a_{2^k+i}(f)=f\Bigl(\frac{2i-1}{2^{k+1}}\Bigr)-\frac12\Bigl(f\Bigl(\frac{i-1}{2^k}\Bigr)+f\Bigl(\frac{i}{2^k}\Bigr)\Bigr).
   $$

   再令 $P_nf(x)=\sum_{j=0}^na_j(f)\varphi_j(x)$（$n=0,1,2,\dots$）。求证：
   (i) 画出 $\varphi_n$（$n=0,\dots,8$）的图像；
   (ii) $m\ne n$ 时 $a_m(\varphi_n)=0$，而 $a_n(\varphi_n)=1$；
   (iii) 若 $n>m$ 则 $P_n\varphi_m=\varphi_m$、$P_m\varphi_n=0$；
   (iv) $P_mP_n=P_{\min\{m,n\}}$；
   (v) $P_0f(0)=f(0)$，且当 $x=0,1$ 或 $x=\frac{2i-1}{2^{r+1}}$（$2^r+i<n$）时 $P_nf(x)=f(x)$；
   (vi) $P_nf$ 是分段线性函数，且在节点 $(0,f(0))$、$(1,f(1))$、$\bigl(\frac{2i-1}{2^{r+1}},f\bigl(\frac{2i-1}{2^{r+1}}\bigr)\bigr)$（$2^r+i<n$）上与 $f$ 相等；
   (vii) 由 $f$ 在 $[0,1]$ 上的一致连续性证明 $P_nf\to f$ 于 $\lVert\cdot\rVert_\infty$，从而 $f(x)=\sum_{n=0}^\infty a_n(f)\varphi_n(x)$；
   (viii) 反之，若 $\sum_{n=0}^\infty b_n\varphi_n$ 在 $C[0,1]$ 中一致收敛，则 $b_n=a_n(f)$ 对一切 $n$ 成立。

> **补充（练习答案）**
>
> **1、2、3 是同一件事。** 这三问只需说明「唯一性」，收敛性由范数的定义直接得到。
>
> - **1.** 设 $x=(\xi_n)\in c_0$、$\varepsilon>0$，取 $N$ 使 $|\xi_n|<\varepsilon$（$n>N$）。令 $\alpha_n=\xi_n$，则 $x-\sum_{i=1}^n\alpha_ie_i$ 只剩尾巴，其 $\lVert\cdot\rVert_\infty<\varepsilon$。唯一性：若 $x=\sum\alpha_ie_i$ 且 $\alpha_n\ne0$，则第 $n$ 个坐标恒为 $\alpha_n\ne0$，与极限在第 $n$ 个坐标等于 $0$ 矛盾（坐标泛函连续）。
> - **2.** 关键是 $c=c_0\oplus\mathbb Ke$，且 $\lVert\cdot\rVert_\infty$ 下 $\alpha_n=\xi_n-\xi_0$（$\xi_0=\lim\xi_n$）、$e$ 的系数 $=\xi_0$，故表示唯一。收敛性由 $\xi_n\to\xi_0$ 得到。
> - **3.** 完全同 1；尾巴由 $\sum_{k>N}|\xi_k|^p\to0$ 控制。
> - **4.** **没有。** 由定理 2.3.4，有 Schauder 基 $\Rightarrow$ 可分，而 $\ell^\infty$ 不可分（例 2.3.2(f)），故无基。
> - **5.** 这就是 **Schauder 系**：把 $[0,1]$ 按二分点不断加密，$\varphi_{2^k+i}$ 是第 $i$ 个 $k$-级子区间 $\bigl(\frac{i-1}{2^k},\frac{i}{2^k}\bigr)$ 上的**帽函数**（峰高 $1$、峰值点 $\frac{2i-1}{2^{k+1}}$），$a_n(f)$ 度量 $f$ 在该峰值点相对两端点连线的偏差。逐条：
>
>   　**(i)** $\varphi_0\equiv1$、$\varphi_1(x)=x$；$\varphi_2,\varphi_3$ 是 $[0,\frac12]$、$[\frac12,1]$ 上的三角帽；$\varphi_4,\dots,\varphi_7$ 是四个四分之一区间上的三角帽；$\varphi_8$ 起是八分之一区间上的。
>
>   　**(ii)** 若 $m$ 与 $n$ 同级而 $m\ne n$：$\varphi_n$ 在该级的**所有**节点上取值 $0$（只有自己的峰值点非零），故 $a_m(\varphi_n)=0$。若 $\mathrm{lev}(m)>\mathrm{lev}(n)$：$\varphi_n$ 在更细的节点上仍为 $0$，同样得 $0$。若 $\mathrm{lev}(m)<\mathrm{lev}(n)$：$\varphi_n$ 的支撑区间 $\bigl(\frac{i-1}{2^k},\frac{i}{2^k}\bigr)$（$k=\mathrm{lev}(n)>m$）严格落在 $m$ 所在区间的**某一半内部**，而帽函数在支撑区间的两个端点取 $0$，故 $m$ 的三个节点（两端与中点）处 $\varphi_n$ 全为 $0$，$a_m(\varphi_n)=0-\frac{0+0}{2}=0$。（这里**不能**说「$\varphi_n$ 在 $m$ 区间上线性」—— 它是个帽子，中间有折点；真正用到的是「三点都取零」。）而 $a_n(\varphi_n)=1-\frac{0+0}{2}=1$（以及 $a_0(\varphi_0)=1$、$a_1(\varphi_1)=1-0=1$）。
>
>   　**(iii)** 由 (ii)，$P_n\varphi_m=\sum_{j\le n}a_j(\varphi_m)\varphi_j$ 只在 $j=m\le n$ 时留下 $\varphi_m$；$n>m$ 时故 $P_n\varphi_m=\varphi_m$，$m<n$ 时 $P_m\varphi_n=0$。
>
>   　**(iv)** 两边都作用在 $\varphi_i$ 上：$i\le\min$ 时都给 $\varphi_i$，$i>\min$ 时都由 (iii) 给 $0$；$\{\varphi_i\}$ 线性无关（各帽函数支撑区间互不包含、按层级严格缩小），故两算子相等。
>
>   　**(v)(vi)** 归纳可见 $P_nf$ 在每个 $k$-级子区间上是线性的，且在 $0,1$ 与已纳入的各峰值点上与 $f$ 相等 —— 它就是 $f$ 在「节点集」上的**分段线性插值**。理由：$a_n(f)\varphi_n$ 只在第 $n$ 个帽的支撑上非零，且在该支撑的两个端点取值 $0$、在峰值点取值 $a_n(f)$，正好把 $f$ 在该点的偏差补上。
>
>   　**(vii)** 设 $k$ 级子区间的长度为 $2^{-k}$。在该区间上 $f$ 与它的线性插值之差不超过 $2\omega_f(2^{-k})$（$\omega_f$ 是连续模），故 $\lVert P_nf-f\rVert_\infty\le2\omega_f\bigl(2^{-k(n)}\bigr)\to0$。于是 $f=\sum_na_n(f)\varphi_n$（按 $\lVert\cdot\rVert_\infty$ 收敛）。
>
>   　**(viii)** 若 $\sum b_n\varphi_n$ 一致收敛到 $f$，则每个 $a_m$ 是有限个点取值（在 $m\ge2$ 时是三点取值）的线性组合，因而在 $(C[0,1],\lVert\cdot\rVert_\infty)$ 上连续，可以逐项作用：$b_m=\sum_nb_na_m(\varphi_n)=a_m(f)$。故系数唯一。
>
>   　这条练习说明**可分性（例 2.3.2(c)）可以升级成「存在 Schauder 基」**，比抽象的稠密性结论强得多。

---

## 订正

1. **例 2.3.2(a)（$c$ 可分）的误差估计差一个因子。** 原文对 $n=0,1,\dots,N$ 取 $|\xi_n-r_n|<\varepsilon$，随后断言 $\lVert x-y\rVert_\infty<\max\{\varepsilon,\ \sup_{n>N}(|\xi_n-\xi_0|+|\xi_0-r_0|)\}=\varepsilon$ 不成立：$|\xi_n-\xi_0|<\varepsilon/2$ 与 $|\xi_0-r_0|<\varepsilon$ 相加只能给出 $\sup_{n>N}|\xi_n-r_0|<3\varepsilon/2$。要凑出 $\lVert x-y\rVert_\infty<\varepsilon$，须取 $|\xi_n-r_n|<\varepsilon/2$（本文正文已按此写）。
2. **定理 2.3.4 中「$r_1,\dots,r_n\in\mathbb Q$」只在实空间成立。** $E$ 是复 Banach 空间时，可数稠密集应取有理**复**数的有限组合，即 $r_i\in\mathbb Q+i\mathbb Q$。本文已在证明中补出这一步。
3. **例 2.3.2(b) 用「$c_0$ 是 $c$ 的子空间故可分」时缺一句理由。** 可分性确实被（任意）子集继承，但这一步不是「显然」：可数稠密集的点一般不在子空间里，须重新挑选。本文在 2.3.2(b) 后补了完整论证，并同时点明「不可分性不被子空间继承」这一反向事实。
4. **定义 2.3.3 的「显然 Schauder 基必线性无关」略去了唯一性的用法。** 本文补出了完整的三行论证。

---

## 下一步

2.4 节转向 Banach 空间中的**完备性与紧性**：级数的绝对收敛刻画完备性、任意赋范空间都有完备化、Riesz 引理与「无限维单位球不紧」、以及 Banach 不动点定理。第 3 章起进入 **Hilbert 空间**：给赋范空间再装一个内积。

---

系列：[[01-预备·Zorn引理与向量空间]] · [[02-度量空间与拓扑空间]] · [[03-拓扑向量空间]] · [[04-赋范空间]] · [[05-Banach空间]] · [[07-Banach空间·完备化与紧性]] · [[08-Hilbert空间·基本性质]] · [[09-Hilbert空间·最佳逼近]] · [[10-Hilbert空间·正交分解]] · [[11-Hilbert空间·正交基]] · [[12-连续线性算子·有界性与算子范数]] · [[13-连续线性算子·连续线性泛函与弱拓扑]] · [[14-连续线性算子·对偶空间的例子]] · [[15-连续线性算子·伴随算子]] · [[16-连续线性算子·投影算子]]
