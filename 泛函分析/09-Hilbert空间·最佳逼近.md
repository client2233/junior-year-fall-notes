# 09 Hilbert 空间 · 最佳逼近

> 系列：[[01-预备·Zorn引理与向量空间]] · [[02-度量空间与拓扑空间]] · [[03-拓扑向量空间]] · [[04-赋范空间]] · [[05-Banach空间]] · [[06-可分Banach空间]] · [[07-Banach空间·完备化与紧性]] · [[08-Hilbert空间·基本性质]] · [[10-Hilbert空间·正交分解]] · [[11-Hilbert空间·正交基]] · [[12-连续线性算子·有界性与算子范数]] · [[13-连续线性算子·连续线性泛函与弱拓扑]] · [[14-连续线性算子·对偶空间的例子]] · [[15-连续线性算子·伴随算子]] · [[16-连续线性算子·投影算子]]
>
> 3.1 建立的是内积「**给出什么**」（范数、夹角、平行四边形律）。3.2 开始「**用**它」：给定 $x\in H$ 与非空闭凸集 $B\subseteq H$，问 $x$ 到 $B$ 的**最短距离是否被某个点取到**。在 Hilbert 空间里答案是肯定的，且取到的点唯一 —— 这就是**投影**的来源，也是 3.3、3.4 乃至第 4 章一切投影算子的地基。
>
> 三件事：
>
> 1. **存在唯一**（定理 3.2.1）：靠平行四边形律 + 完备性，用一条**极小化序列**把最优点造出来。
> 2. **变分刻画**（定理 3.2.2）：把「$\tilde x$ 是最佳逼近元」翻译成一个**内积不等式**，于是「最优」成了可计算的条件。
> 3. **投影映射非扩张**（推论 3.2.6）：最佳逼近点连续依赖 $x$，且不放大距离。

## 3.2 最佳逼近

### 3.2.1 距离与最佳逼近元

设 $X$ 是度量空间，$\varnothing\ne B\subseteq X$，$x\in X\setminus B$。$x$ 到 $B$ 的**距离**定义为

$$
d(x,B)=\inf_{y\in B}d(x,y).
$$

若存在 $\tilde x\in B$ 使 $d(x,\tilde x)=d(x,B)$，则称 $\tilde x$ 是 $x$ 从 $B$ 的**最佳逼近元**（best approximate element）。

**要点（两条与 $B$ 形状无关的一般事实）**

- $d(x,B)=0$ 当且仅当 $x\in\overline B$。特别地，$B$ 闭时 $d(x,B)=0\iff x\in B$。
- 对固定 $B$，函数 $x\mapsto d(x,B)$ 是 $1$-Lipschitz 的：

$$
\lvert d(x,B)-d(y,B)\rvert\le d(x,y).
$$

  这是 $\inf$ 的一般性质（任取 $b$ 有 $d(x,B)\le d(x,y)+d(y,b)$，对 $b$ 取下确界即得），与 $B$ 的形状无关。

**要点（$\inf$ 与 $\min$ 的差别就是本节全部内容）** 「距离」由 $\inf$ 定义，**总存在**（下确界 axiomatically 存在）；「最佳逼近元」要求这个 $\inf$ 被某个 $B$ 中的点**取到**，而这不总成立。定理 3.2.1 说的就是：在 Hilbert 空间里、对闭凸集，$\inf$ 一定取到、且取到的点唯一。三个条件（Hilbert 性、闭、凸）缺一不可，3.2.5 给出三组反例。

### 3.2.2 存在唯一性（定理 3.2.1）

**定理 3.2.1** 设 $B$ 是 Hilbert 空间 $H$ 的非空闭凸子集，$x\in H\setminus B$。则存在唯一的 $\tilde x\in B$ 使

$$
\lVert x-\tilde x\rVert=d(x,B)=\inf_{y\in B}\lVert x-y\rVert.
$$

**证明（一）存在性：造极小化序列。** 记 $d=d(x,B)$。由 $\inf$ 的定义，对每个 $n\in\mathbb N$ 存在 $y_n\in B$ 使

$$
d\le\lVert x-y_n\rVert<d+\frac1n.
$$

断言 $(y_n)_{n\in\mathbb N}$ 是 Cauchy 列。对 $u=x-y_n$、$v=x-y_m$ 用平行四边形律（定理 3.1.10，即 $\lVert u+v\rVert^2+\lVert u-v\rVert^2=2\lVert u\rVert^2+2\lVert v\rVert^2$），注意

$$
u+v=2x-y_n-y_m,\qquad u-v=y_m-y_n,
$$

得

$$
2\lVert x-y_n\rVert^2+2\lVert x-y_m\rVert^2=\lVert 2x-y_n-y_m\rVert^2+\lVert y_n-y_m\rVert^2.
$$

由 $B$ 凸，$\dfrac{y_n+y_m}{2}\in B$，于是

$$
\lVert 2x-y_n-y_m\rVert^2=4\Bigl\lVert x-\frac{y_n+y_m}{2}\Bigr\rVert^2\ge4d^2.
$$

代入并整理：

$$
\lVert y_n-y_m\rVert^2\le2\Bigl(d+\frac1n\Bigr)^2+2\Bigl(d+\frac1m\Bigr)^2-4d^2=\frac{4d}{n}+\frac{4d}{m}+\frac{2}{n^2}+\frac{2}{m^2}\longrightarrow0\qquad(m,n\to\infty).
$$

故 $(y_n)$ 是 Cauchy 列。由 $H$ 完备、$B$ 闭，存在 $\tilde x\in B$ 使 $y_n\to\tilde x$。再由范数的连续性，

$$
\lVert x-\tilde x\rVert=\lim_{n\to\infty}\lVert x-y_n\rVert=d,
$$

即 $\tilde x$ 是 $x$ 从 $B$ 的最佳逼近元。

**（二）唯一性。** 设 $\tilde x_1\in B$ 也满足 $\lVert x-\tilde x_1\rVert=d$。对 $u=x-\tilde x$、$v=x-\tilde x_1$ 再用一次平行四边形律（此时 $u+v=2x-\tilde x-\tilde x_1$、$u-v=\tilde x_1-\tilde x$）：

$$
4d^2=2\lVert x-\tilde x\rVert^2+2\lVert x-\tilde x_1\rVert^2=\lVert\tilde x-\tilde x_1\rVert^2+4\Bigl\lVert x-\frac{\tilde x+\tilde x_1}{2}\Bigr\rVert^2\ge\lVert\tilde x-\tilde x_1\rVert^2+4d^2.
$$

末步用了 $\dfrac{\tilde x+\tilde x_1}{2}\in B$（凸性）。于是 $\lVert\tilde x-\tilde x_1\rVert^2\le0$，即 $\tilde x=\tilde x_1$。$\square$

**要点（平行四边形律在这里究竟做了什么）** 一般赋范空间里，极小化序列**不必**是 Cauchy 列 —— 「$B$ 的两端各自逼近 $x$」完全可能，此时点列来回跳动而不收敛。平行四边形律恰好把这件事堵死，因为它把 $\lVert y_n-y_m\rVert^2$ 解了出来：

$$
\lVert y_n-y_m\rVert^2=2\lVert x-y_n\rVert^2+2\lVert x-y_m\rVert^2-4\Bigl\lVert x-\frac{y_n+y_m}{2}\Bigr\rVert^2.
$$

右端最后一项被 $4d^2$ 从下方压住 —— **这是凸性被用到的唯一位置**；前两项趋于 $2d^2$。没有这条恒等式，右端消不掉。这是「内积范数」比一般范数好用的第一处实质体现（第二处是 3.3 的正交分解）。

**补充（这条定理是一般定理的特例）** 定理 3.2.1 真正需要的只是「$H$ 是**一致凸**（uniformly convex）的 Banach 空间」：对任意 $\varepsilon>0$ 存在 $\delta>0$，使 $\lVert u\rVert,\lVert v\rVert\le1$ 且 $\lVert u-v\rVert\ge\varepsilon$ 时可推出 $\lVert(u+v)/2\rVert\le1-\delta$。Hilbert 空间由平行四边形律直接给出一致凸性（取 $\delta=1-\sqrt{1-\varepsilon^2/4}$ 即可）。一致凸蕴含**严格凸**，严格凸给出**唯一性**；再配上自反性就给出**存在性**（Milman–Pettis 定理）。所以本节这条定理是「Hilbert 空间是三类好空间（一致凸 / 自反 / 内积）中最具体的一个」这句话的起点。

**补充（假设里的「闭」在定理 3.2.2 里其实是多余的）** 定理 3.2.2 的等价性只用凸性；闭性只在定理 3.2.1 的**存在性**里被用到（保证 Cauchy 列的极限留在 $B$ 内）。课件与教材都把「闭」一起写进 3.2.2 的假设，读的时候知道它是冗余的即可。

### 3.2.3 变分刻画（定理 3.2.2）

**定理 3.2.2** 设 $B$ 是 Hilbert 空间 $H$ 中的闭凸集，$\tilde x\in B$、$x\in H\setminus B$。以下三条等价：

(1) $\lVert x-\tilde x\rVert=d(x,B)$；

(2) 对一切 $b\in B$，$\operatorname{Re}\langle x-\tilde x,\,b-\tilde x\rangle\le0$；

(3) 对一切 $b\in B$，$\operatorname{Re}\langle x-b,\,\tilde x-b\rangle\ge0$。

**证明。**

$(1)\Longrightarrow(2)$。对任意 $b\in B$ 与 $0<\lambda<1$，由凸性

$$
(1-\lambda)\tilde x+\lambda b=\tilde x+\lambda(b-\tilde x)\in B.
$$

于是

$$
\lVert x-\tilde x\rVert^2=d(x,B)^2\le\bigl\lVert x-\tilde x-\lambda(b-\tilde x)\bigr\rVert^2=\lVert x-\tilde x\rVert^2-2\lambda\operatorname{Re}\langle x-\tilde x,\,b-\tilde x\rangle+\lambda^2\lVert b-\tilde x\rVert^2.
$$

消去首项、两边除以 $\lambda>0$：

$$
2\operatorname{Re}\langle x-\tilde x,\,b-\tilde x\rangle\le\lambda\lVert b-\tilde x\rVert^2.
$$

令 $\lambda\to0^+$ 即得 (2)。

$(2)\Longrightarrow(3)$。把 $x-b$ 拆成 $(x-\tilde x)+(\tilde x-b)$，对 $b\in B$ 用 (2)：

$$
\operatorname{Re}\bigl\langle(x-b)+(b-\tilde x),\,b-\tilde x\bigr\rangle\le0,
$$

即

$$
\operatorname{Re}\langle x-b,\,b-\tilde x\rangle+\lVert b-\tilde x\rVert^2\le0,
$$

故 $\operatorname{Re}\langle x-b,\,b-\tilde x\rangle\le0$，等价于 $\operatorname{Re}\langle x-b,\,\tilde x-b\rangle\ge0$，即 (3)。

$(3)\Longrightarrow(2)$。对 $b\in B$、$0<\lambda<1$，令 $b_\lambda=\tilde x+\lambda(b-\tilde x)\in B$，对 $b_\lambda$ 用 (3)：

$$
\operatorname{Re}\bigl\langle x-b_\lambda,\,\tilde x-b_\lambda\bigr\rangle=\operatorname{Re}\bigl\langle x-\tilde x-\lambda(b-\tilde x),\,-\lambda(b-\tilde x)\bigr\rangle\ge0,
$$

即

$$
-\lambda\operatorname{Re}\langle x-\tilde x,\,b-\tilde x\rangle+\lambda^2\lVert b-\tilde x\rVert^2\ge0.
$$

除以 $\lambda$、令 $\lambda\to0^+$，得 (2)。

$(2)\Longrightarrow(1)$。对任意 $b\in B$，

$$
\lVert x-b\rVert^2=\bigl\lVert(x-\tilde x)+(\tilde x-b)\bigr\rVert^2=\lVert x-\tilde x\rVert^2+2\operatorname{Re}\langle x-\tilde x,\,\tilde x-b\rangle+\lVert\tilde x-b\rVert^2\ge\lVert x-\tilde x\rVert^2,
$$

其中用了 $\operatorname{Re}\langle x-\tilde x,\,\tilde x-b\rangle=-\operatorname{Re}\langle x-\tilde x,\,b-\tilde x\rangle\ge0$。对 $b$ 取下确界即得 $\lVert x-\tilde x\rVert=d(x,B)$。$\square$

**要点（几何意义：支撑超平面）** 条件 (2) 说的是：**$B$ 整个落在过 $\tilde x$、以 $x-\tilde x$ 为外法向的闭半空间内**。等价地说，以 $x$ 为心、$d$ 为半径的闭球与 $B$ 在 $\tilde x$ 处「相切」（只接触不穿入）。这就把「求最短距离」变成「找一个支撑超平面的切点」。条件 (3) 是同一件事的另一种写法：在三角形 $b\tilde x x$ 的顶点 $b$ 处，角 $\angle xb\tilde x$ 是锐角；两式相减恰好余下 $\lVert b-\tilde x\rVert^2\ge0$。

**要点（$\operatorname{Re}$ 的作用）** 复空间里 $\langle x-\tilde x,\,b-\tilde x\rangle$ 是复数，「$\le0$」无从谈起，必须取实部；实空间（$K=\mathbb R$）里 $\operatorname{Re}$ 可以整体删去（注 3.2.3）。

**补充（(2) 就是一阶最优性条件）** 令 $\varphi(\lambda)=\lVert x-\tilde x-\lambda(b-\tilde x)\rVert^2$，则 $\tilde x$ 是极小点当且仅当右导数 $\varphi'(0^+)\ge0$，而 $\varphi'(0^+)=-2\operatorname{Re}\langle x-\tilde x,\,b-\tilde x\rangle$。所以 (2) 就是无约束极值里 $\phi'(0^+)\ge0$ 那个条件的精确版本 —— 区别只在于这里的「可行方向」$b-\tilde x$ 是**整个凸集 $B$ 给出的锥**，不是切空间。

### 3.2.4 投影映射（推论 3.2.4、定义 3.2.5、推论 3.2.6）

**推论 3.2.4** 设 $B$ 是 Hilbert 空间 $H$ 的非空闭凸子集，$x,y\in H$，$\tilde x,\tilde y$ 分别是 $x,y$ 从 $B$ 的最佳逼近元，则

$$
\lVert\tilde x-\tilde y\rVert\le\lVert x-y\rVert.
$$

**证明。** 对任意 $b\in B$，由定理 3.2.2 的 (2) 有

$$
\operatorname{Re}\langle x-\tilde x,\,b-\tilde x\rangle\le0,\qquad\operatorname{Re}\langle y-\tilde y,\,b-\tilde y\rangle\le0.
$$

分别取 $b=\tilde y$ 与 $b=\tilde x$，相加：

$$
\operatorname{Re}\langle x-\tilde x,\,\tilde y-\tilde x\rangle+\operatorname{Re}\langle y-\tilde y,\,\tilde x-\tilde y\rangle\le0.
$$

两式的被加项合并成一项：注意 $\tilde y-\tilde x=-(\tilde x-\tilde y)$，故

$$
\langle x-\tilde x,\,\tilde y-\tilde x\rangle+\langle y-\tilde y,\,\tilde x-\tilde y\rangle=\bigl\langle-(x-\tilde x)+(y-\tilde y),\,\tilde x-\tilde y\bigr\rangle=\bigl\langle(y-x)+(\tilde x-\tilde y),\,\tilde x-\tilde y\bigr\rangle,
$$

取实部后不等式成为

$$
\operatorname{Re}\langle y-x,\,\tilde x-\tilde y\rangle+\lVert\tilde x-\tilde y\rVert^2\le0,
$$

即

$$
\lVert\tilde x-\tilde y\rVert^2\le\operatorname{Re}\langle x-y,\,\tilde x-\tilde y\rangle .
$$

于是由 Cauchy–Schwarz 不等式（定理 3.1.3）

$$
\lVert\tilde x-\tilde y\rVert^2\le\bigl\lvert\langle x-y,\,\tilde x-\tilde y\rangle\bigr\rvert\le\lVert x-y\rVert\,\lVert\tilde x-\tilde y\rVert .
$$

若 $\tilde x\ne\tilde y$，两边约去 $\lVert\tilde x-\tilde y\rVert>0$ 即得结论；若 $\tilde x=\tilde y$ 则结论显然。$\square$

**定义 3.2.5** 设 $B$ 是 Hilbert 空间 $H$ 的非空闭凸子集。定义**投影映射** $P_B:H\to B$ 为

$$
P_B(x)=\tilde x,
$$

其中 $\tilde x$ 是 $x$ 从 $B$ 的唯一最佳逼近元。

**推论 3.2.6** 投影映射 $P_B$ 连续；事实上它是 $1$-Lipschitz 的：

$$
\lVert P_Bx-P_By\rVert\le\lVert x-y\rVert,\qquad x,y\in H.
$$

**要点（$P_B$ 什么时候是线性的）** 一般情形下 $\tilde x$ 依赖 $x$ 是**非线性**的。判据（练习 4）：$P_B$ 线性 $\iff$ $B$ 是 $H$ 的线性子空间。当 $B=M$ 是闭子空间时，$P_M$ 不仅线性，还**幂等**（$P_M^2=P_M$）、**自伴**（$\langle P_Mx,y\rangle=\langle x,P_My\rangle$），并且 $\lVert P_M\rVert\le1$。这一族性质正是第 4 章 4.5 节「Hilbert 空间上的投影」的主题。3.3 先把「$B=M$ 是闭子空间」这一特殊情形单独做成定理（定理 3.3.5、3.3.7）。

### 3.2.5 三个反例：Hilbert 性、闭、凸缺一不可

**（一）闭性不能省（练习 2）** 取 $H=\ell^2$，

$$
B=\Bigl\{(\xi_n)_{n\in\mathbb N}\in\ell^2:\ \xi_n\ne0\ \text{只有有限多个},\ \sum_{n=1}^\infty\xi_n=1\Bigr\}.
$$

$B$ 凸（两个有限支集序列的凸组合仍是有限支集，求和仍为 $1$），但 $B$ **不闭**：$y_n=\bigl(\underbrace{\tfrac1n,\dots,\tfrac1n}_{n},0,0,\dots\bigr)$ 都在 $B$ 中，而

$$
\lVert y_n\rVert_2=\frac1{\sqrt n}\longrightarrow0.
$$

于是 $y_n\to0\notin B$。此时 $d(0,B)=0$ —— 因为 $\lVert y_n\rVert\to0$ —— 但 $0\notin B$，所以 $0$ 在 $B$ 中**没有**最佳逼近元：要 $\lVert 0-y\rVert=0$ 只能 $y=0$。**这正是定理 3.2.1 里「$B$ 闭」的用处：没有闭性，极限点会跑出 $B$。**

**（二）凸性不能省（补充反例）** 取 $H=\mathbb R^2$（Euclidean 范数）、$B=\{(x_1,x_2):x_1^2+x_2^2=1\}$（闭但不凸）、$x=0$。则 $d(0,B)=1$，但 $B$ 中每一点到 $0$ 的距离都等于 $1$ —— 最佳逼近元**处处都是**，唯一性彻底失败。可见「凸」是唯一性的直接来源。

**（三）Hilbert 性不能省（练习 3）**

**(a) 存在性会失败。** 取 $X=C[-1,1]$ 配 $\sup$ 范数（Banach 空间，但不是内积空间 —— $\sup$ 范数不满足平行四边形律），

$$
M=\Bigl\{f\in C[-1,1]:\ \int_{-1}^0f(t)\,dt=\int_0^1f(t)\,dt=0\Bigr\}
$$

是两个连续线性泛函的公共零点，故是闭线性子空间（从而闭凸）。取 $g\in C[-1,1]$ 满足

$$
\int_{-1}^0g(t)\,dt=1,\qquad \int_0^1g(t)\,dt=-1.
$$

对任意 $f\in M$，$g-f$ 的两个积分分别仍为 $1$ 与 $-1$，而区间长度都是 $1$，故

$$
\lVert g-f\rVert_\infty\ge\Bigl\lvert\int_{-1}^0(g-f)\Bigr\rvert=1,\qquad \lVert g-f\rVert_\infty\ge\Bigl\lvert\int_0^1(g-f)\Bigr\rvert=1,
$$

即 $d(g,M)\ge1$。

另一方面 $d(g,M)\le1$。对 $\eta>0$ 取分段线性连续函数 $h_\eta$：在 $[-1,-\eta]$ 上恒等于 $A$，在 $[\eta,1]$ 上恒等于 $-A$，在 $[-\eta,\eta]$ 上线性过渡（于是 $h_\eta(0)=\frac{A+(-A)}{2}=0$），其中

$$
A=\frac{1}{1-\eta/2}.
$$

注意 $\lVert h_\eta\rVert_\infty=A$。逐段算积分（$[-\eta,0]$ 与 $[0,\eta]$ 上分别是「$A$ 到 $0$」「$0$ 到 $-A$」的线性函数，各贡献长度乘平均高度）：

$$
\int_{-1}^0h_\eta=(1-\eta)A+\frac{\eta}{2}\,(A+0)=A\Bigl(1-\frac\eta2\Bigr)=1,
$$

$$
\int_0^1h_\eta=\frac{\eta}{2}\,\bigl(0+(-A)\bigr)+(1-\eta)(-A)=-A\Bigl(1-\frac\eta2\Bigr)=-1 .
$$

于是 $f_\eta:=g-h_\eta$ 满足 $\int_{-1}^0f_\eta=1-1=0$、$\int_0^1f_\eta=-1-(-1)=0$，即 $f_\eta\in M$，而

$$
\lVert g-f_\eta\rVert_\infty=\lVert h_\eta\rVert_\infty=A=\frac{1}{1-\eta/2}\longrightarrow1\qquad(\eta\to0^+).
$$

故 $d(g,M)=1$，但**取不到**：若 $f\in M$ 使 $\lVert g-f\rVert_\infty=1$，则 $\lvert g-f\rvert\le1$ 且 $\int_{-1}^0(g-f)=1$，由长度 $1$ 与连续性迫使 $g-f\equiv1$ 于 $[-1,0]$；同理 $g-f\equiv-1$ 于 $[0,1]$ —— 在 $t=0$ 处矛盾。

**(b) 唯一性会失败。** 取 $X=\mathbb R^2$ 配 $\lVert(\xi_1,\xi_2)\rVert_\infty=\max\{\lvert\xi_1\rvert,\lvert\xi_2\rvert\}$（不是内积范数），$B=\{(\xi_1,\xi_2)\in X:\xi_1\ge1\}$（闭凸）。原点到 $B$ 的距离为 $1$，而线段

$$
\{(\xi_1,\xi_2)\in X:\ \xi_1=1,\ \lvert\xi_2\rvert\le1\}
$$

上**每一点**到原点的 $\infty$-距离都等于 $1$。最佳逼近元存在但不唯一。

**小结** 三个反例对应定理 3.2.1 的三个假设，逻辑很干净：

| 去掉的假设 | 反例 | 失败的是 |
|---|---|---|
| 闭 | $\ell^2$ 中「有限支集且求和为 $1$」的集合 | 存在性 |
| 凸 | $\mathbb R^2$ 中的单位圆周 | 唯一性 |
| Hilbert（内积范数） | (a) $C[-1,1]$ 的 $\sup$ 范数；(b) $\mathbb R^2$ 的 $\infty$-范数 | (a) 存在性；(b) 唯一性 |

### 3.2.6 练习

教材 3.2 共四题：(1)(4) 的解答写在下面，(2)(3) 是定理 3.2.1 三个假设的反例**验证题**，完整计算已在 §3.2.5 给出，此处只标出对应位置与一句话结论。

**练习 1** 设 $K_1,K_2$ 是 Hilbert 空间 $H$ 的非空闭凸子集，$K_1\subseteq K_2$。对任意 $x\in H$，设 $y_1,y_2$ 分别是 $x$ 在 $K_1,K_2$ 中的最近点，证明

$$
\lVert y_1-y_2\rVert^2\le2\bigl(d(x,K_1)^2-d(x,K_2)^2\bigr).
$$

**解。** 记 $d_i=d(x,K_i)=\lVert x-y_i\rVert$（$i=1,2$）。由 $K_1\subseteq K_2$ 得 $d_1\ge d_2$，故右端非负。对 $u=x-y_1$、$v=x-y_2$ 用平行四边形律：

$$
\lVert y_1-y_2\rVert^2=2\lVert x-y_1\rVert^2+2\lVert x-y_2\rVert^2-4\Bigl\lVert x-\frac{y_1+y_2}{2}\Bigr\rVert^2.
$$

关键一步：$y_1\in K_1\subseteq K_2$、$y_2\in K_2$，而 $K_2$ 凸，所以 $\dfrac{y_1+y_2}{2}\in K_2$，于是

$$
\Bigl\lVert x-\frac{y_1+y_2}{2}\Bigr\rVert\ge d_2.
$$

代入即得

$$
\lVert y_1-y_2\rVert^2\le2d_1^2+2d_2^2-4d_2^2=2(d_1^2-d_2^2).\qquad\square
$$

**要点** 这条不等式说的是：**$K$ 变小 $d(x,K)$ 变大，而「最优点的移动量」被距离的增量平方根控制**。它是「最佳逼近对约束集的连续依赖」的定量版本。

**练习 2（闭性）** 取 $H=\ell^2$，$B$ 为「只有有限多个非零坐标、且 $\sum_{n\ge1}\xi_n=1$」的序列全体。证明 $B$ 凸、不闭，且 $0$ 在 $B$ 中没有最佳逼近元。**答。** $y_n=(\frac1n,\dots,\frac1n,0,\dots)\in B$ 而 $\lVert y_n\rVert_2=1/\sqrt n\to0$，故 $0\in\overline B\setminus B$，$d(0,B)=0$；但 $\lVert y\rVert_2=0$ 只有 $y=0\notin B$，所以距离取不到。完整论证见 §3.2.5「（一）闭性不能省」。

**练习 3（Hilbert 性）** $(a)$ 在 $(C[-1,1],\lVert\cdot\rVert_\infty)$ 中取 $M=\{f:\int_{-1}^0f=\int_0^1f=0\}$、$g$ 满足两个积分分别为 $1$ 与 $-1$。证明 $d(g,M)=1$ 但取不到。$(b)$ 在 $(\mathbb R^2,\lVert\cdot\rVert_\infty)$ 中取 $B=\{\xi_1\ge1\}$，证明原点在 $B$ 中的最佳逼近元不唯一。**答。** $(a)$ 对 $f\in M$ 有 $\lVert g-f\rVert_\infty\ge\lvert\int_{-1}^0(g-f)\rvert=1$，而用斜坡函数 $h_\eta$ 从上方逼近可得 $d(g,M)=1$；若取到，则 $g-f$ 在两半区间上分别恒为 $\pm1$，在 $t=0$ 处矛盾。$(b)$ $d(0,B)=1$，整条线段 $\{\xi_1=1,\lvert\xi_2\rvert\le1\}$ 上每点都达到 $1$。完整论证见 §3.2.5「（三）Hilbert 性不能省」。

**练习 4** 证明：非空闭凸子集 $B\subseteq H$ 的投影映射 $P_B$ 线性，当且仅当 $B$ 是 $H$ 的线性子空间。

**解。**

$(\Longleftarrow)$ 设 $B=M$ 是闭线性子空间。由 3.3 的定理 3.3.5，$P_Mx$ 是唯一的 $y\in M$ 使 $(x-y)\perp M$；这条刻画对 $x$ 线性（唯一性由正交分解保证），故 $P_M$ 线性。

$(\Longrightarrow)$ 设 $P_B$ 线性。$P_B$ 是满射（$b\in B$ 时 $P_Bb=b$），故 $B=P_B(H)$ 是线性映射的像，即线性子空间。$\square$

**要点** 这条练习解释了为什么 3.3 要把「$B$ 是闭子空间」单独拿出来做：只有这种 $B$，投影才是一步线性代数（解方程组），而不是一个变分问题。

---

## 订正

1. **定理 3.2.1 证明里有一处 $\le$ 应为 $=$。** 原文写

   $\lVert x-\tilde x\rVert^2\le\lVert x-(\tilde x+\lambda(b-\tilde x))\rVert^2\le\lVert x-\tilde x\rVert^2-2\operatorname{Re}\langle x-\tilde x,\lambda(b-\tilde x)\rangle+\lvert\lambda\rvert^2\lVert b-\tilde x\rVert^2$。

   第二个不等号是错的：中间那项按内积展开是**恒等式**。本文已改为 $=$。（第一个 $\le$ 是对的，它来自 $\lVert x-\tilde x\rVert=d(x,B)$。）
2. **推论 3.2.6 把 $1$-Lipschitz 说成「contraction（压缩）」不准确。** 压缩通常要求 Lipschitz 常数 $<1$（或 $\le k<1$，如定理 2.4.9 的 Banach 不动点定理所要求）。$L=1$ 的映射通行名称是**非扩张**（nonexpansive）。本文用「非扩张 / $1$-Lipschitz」。
3. **定理 3.2.1 与 3.2.2 的证明把定理 3.1.10 称作「parallelism law」。** 该定理的标题是 `Parallelogram Law`，正文中却出现两种叫法（parallelism / parallelogram）。同一律两名，本文统一为「平行四边形律」。
4. **定理 3.2.2 的假设中「闭」是多余的**（见 3.2.2 补充），原文未说明。

---

## 下一步

3.3 把本节的结果用到**闭子空间**上：这时变分不等式 (2) 会升级成等式 $(x-y)\perp M$，于是整个 Hilbert 空间可以劈成 $H=M\oplus M^\perp$ —— 这就是**正交分解定理**，也是第 4 章投影算子、正交基（3.4）与 Fourier 展开的共同地基。

---

系列：[[01-预备·Zorn引理与向量空间]] · [[02-度量空间与拓扑空间]] · [[03-拓扑向量空间]] · [[04-赋范空间]] · [[05-Banach空间]] · [[06-可分Banach空间]] · [[07-Banach空间·完备化与紧性]] · [[08-Hilbert空间·基本性质]] · [[10-Hilbert空间·正交分解]] · [[11-Hilbert空间·正交基]] · [[12-连续线性算子·有界性与算子范数]] · [[13-连续线性算子·连续线性泛函与弱拓扑]] · [[14-连续线性算子·对偶空间的例子]] · [[15-连续线性算子·伴随算子]] · [[16-连续线性算子·投影算子]]
