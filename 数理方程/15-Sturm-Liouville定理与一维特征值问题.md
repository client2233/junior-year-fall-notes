# 15 Sturm–Liouville 定理与一维特征值问题

> 系列：[[14-谱理论·鼓面的数学建模与Laplace特征值]] · [[16-特征值的变分刻画与Poincaré不等式]] · [[12-习题课]]

高维的特征值问题讲起来容易，算起来几乎不可能 —— 除了矩形和圆盘（[[17-矩形与圆盘上的分离变量·Bessel方程]]），一般区域写不出显式解。**一维是唯一能算到底的情形**，而它的答案（三角函数 + 离散谱）正是后面 [[18-热方程·建模、边值条件与分离变量法]] 里 Fourier 级数的来源。

本篇把两讲的同一件事合并：讲稿第六周 §3.3 给了定理的完整陈述，第七周 §1 给了边值条件的分类与**相位法求谱**。两者的编号、措辞都不一致（第六周说「所有特征值 $\lambda\ge0$」，第七周又重复陈述一遍），本篇按知识只留一套，把两边的独有内容并进去。

---

## 15.1 问题与边值条件的分类

令 $\Omega=(0,L)\subseteq\mathbb R$，考虑

$$
\begin{cases}-u''(x)=\lambda u(x), & x\in(0,L),\\[4pt] -\alpha_1u'(0)+\beta_1u(0)=0,\\[4pt] \alpha_2u'(L)+\beta_2u(L)=0,\end{cases}\tag{15.1}
$$

其中

$$
\alpha_i\ge0,\qquad\beta_i\ge0,\qquad\alpha_i+\beta_i\ne0\quad(i=1,2).
$$

**三类边值的判别**

| 条件 | 边值类型 | 说明 |
|---|---|---|
| $\alpha_1=\alpha_2=0$ | Dirichlet | 此时必 $\beta_i\ne0$，即 $u(0)=u(L)=0$ |
| $\beta_1=\beta_2=0$ | Neumann | 此时必 $\alpha_i\ne0$，即 $u'(0)=u'(L)=0$ |
| $\beta_1/\alpha_1>0$ 且 $\beta_2/\alpha_2>0$ | Robin | 比值就是 $\sigma$，见下 |

> **要点（符号从哪来）**　一维区域 $(0,L)$ 上，外法向在左端点指向 $-x$、右端点指向 $+x$，所以 $\partial_\nu u|_{x=0}=-u'(0)$、$\partial_\nu u|_{x=L}=u'(L)$。于是 Robin 条件 $\partial_\nu u+\sigma u=0$ 在两端分别化为 $-u'(0)+\sigma u(0)=0$ 与 $u'(L)+\sigma u(L)=0$，与 (15.1) 逐项对比即得 $\sigma=\beta_1/\alpha_1=\beta_2/\alpha_2$。**两端取同一个 $\sigma$ 只是一种特例**；一般允许两端不同的 $\sigma_1=\beta_1/\alpha_1$、$\sigma_2=\beta_2/\alpha_2$，(15.1) 已经涵盖。

---

## 15.2 Sturm–Liouville 定理

**定理 15.2.1（Sturm–Liouville）**　对问题 (15.1)：

1. 所有特征值 $\lambda\ge0$，且可排成

$$
0\le\lambda_1\le\lambda_2\le\cdots\le\lambda_n\le\cdots\longrightarrow+\infty .
$$

特别地，$\lambda_1=0$ **当且仅当** $\beta_1=\beta_2=0$（即 Neumann 情形）。

2. 存在一列特征函数 $\{u_k\}_{k=1}^\infty\subseteq C^\infty([0,L])$，构成 $L^2((0,L))$ 的完备标准正交基。对任意 $f\in L^2((0,L))$，

$$
f=\sum_{k=1}^\infty a_ku_k(x),\qquad a_k=\int_0^Lf(x)u_k(x)\,dx ,
$$

左右两端在 $L^2((0,L))$ 意义下相等。

3. 若 $f\in C^2([0,L])$ 且 $f$ 满足 (15.1) 的边值条件，则上式右端绝对收敛、且一致收敛到 $f$。

**证明（第 1 条）**　设 $u$ 是 $\lambda$ 对应的特征函数。两边乘 $u$ 并在 $(0,L)$ 上积分，再用分部积分拆成「端点项 + 平方项」：

$$
\lambda\int_0^Lu^2\,dx=\int_0^Lu\,(-u'')\,dx=\int_0^L-(uu')'\,dx+\int_0^L(u')^2\,dx ,
$$

$$
\lambda\int_0^Lu^2\,dx=-\bigl[u(L)u'(L)-u(0)u'(0)\bigr]+\int_0^L(u')^2\,dx .
$$

现在用边值条件把端点项化掉。设 $\alpha_2>0$，由 $\alpha_2u'(L)+\beta_2u(L)=0$ 得 $u'(L)=-\dfrac{\beta_2}{\alpha_2}u(L)$，故 $-u(L)u'(L)=\dfrac{\beta_2}{\alpha_2}u(L)^2\ge0$；同理（$\alpha_1>0$ 时）由 $-\alpha_1u'(0)+\beta_1u(0)=0$ 得 $u'(0)=\dfrac{\beta_1}{\alpha_1}u(0)$，故 $u(0)u'(0)=\dfrac{\beta_1}{\alpha_1}u(0)^2\ge0$。若某个 $\alpha_i=0$，则相应端点的 $u$ 或 $u'$ 本身为零、端点项直接为 $0$，可看作上式的极限情形。代回得

$$
\lambda\int_0^Lu^2\,dx=\frac{\beta_2}{\alpha_2}u(L)^2+\frac{\beta_1}{\alpha_1}u(0)^2+\int_0^L(u')^2\,dx\ge0 .
$$

分母 $\int u^2>0$（$u\not\equiv0$），故 $\lambda\ge0$。

再看 $\lambda=0$。上式给出 $\lambda=0$ 时右端三项全为 $0$，故 $u'\equiv0$、$u\equiv$ 常数。此时端点项 $(\beta_2/\alpha_2)u(L)^2+(\beta_1/\alpha_1)u(0)^2=0$ 要求 $\beta_1=\beta_2=0$（常数非零使 $u(0)=u(L)\ne0$）；反之若 $\beta_1=\beta_2=0$，常函数确实满足两个边值条件，$\lambda=0$ 是特征值。第 1 条得证。$\square$

**补充（第 2、3 条为什么不在这里证）**　它们属于 Sturm–Liouville 理论的标准结论，讲稿直接引用。机制与 [[14-谱理论·鼓面的数学建模与Laplace特征值]] §14.5 的补充完全一样：把 $-\dfrac{d^2}{dx^2}$ 配边值条件后的逆算子写出来（就是用 Green 函数积出来的那个积分算子 $(Tf)(x)=\int_0^LG(x,y)f(y)\,dy$），它是 $L^2(0,L)$ 上的**紧、自伴、正**算子，于是紧自伴算子的谱定理给出：特征值 $\mu_k=1/\lambda_k\searrow0$、每个特征空间有限维、特征向量构成完备正交基。第 3 条的一致收敛还要用 $G$ 的正则性与 $\lambda_k\sim(k\pi/L)^2$ 换来的衰减：$a_k=O(1/\lambda_k)$ 保证 $\sum|a_k|<\infty$，而 $\lambda_ka_k$ 是 $f''$ 的 Fourier 系数、有界。

---

## 15.3 特征值的显式计算（相位法）

定理只保证 $\lambda$ 存在、可数、趋于无穷；要**算出来**得解方程。已证 $\lambda\ge0$，分 $\lambda=0$ 与 $\lambda>0$。

**$\lambda=0$**：$u\equiv$ 常数，由 §15.2 的论证只有 Neumann 时合法，此时 $\lambda_1=0$。

**$\lambda>0$**：通解

$$
u(x)=A\cos(\sqrt\lambda\,x)+B\sin(\sqrt\lambda\,x).
$$

**第一步：代入 $x=0$ 的边值条件。** 由 $u(0)=A$、$u'(0)=B\sqrt\lambda$，

$$
-\alpha_1\sqrt\lambda\,B+\beta_1A=0 .
$$

$\beta_1>0$ 时可取 $A=\alpha_1\sqrt\lambda$、$B=\beta_1$（整体乘常数不影响），得

$$
u(x)=C\bigl[\alpha_1\sqrt\lambda\cos(\sqrt\lambda\,x)+\beta_1\sin(\sqrt\lambda\,x)\bigr]=C\,R_1\sin\bigl(\sqrt\lambda\,x+\theta_1\bigr),
$$

其中

$$
R_1=\sqrt{\alpha_1^2\lambda+\beta_1^2},\qquad
\sin\theta_1=\frac{\alpha_1\sqrt\lambda}{R_1},\quad\cos\theta_1=\frac{\beta_1}{R_1},
$$

即

$$
\theta_1=\begin{cases}\arctan\left(\dfrac{\alpha_1\sqrt\lambda}{\beta_1}\right), & \beta_1>0,\\[8pt] \dfrac{\pi}{2}, & \beta_1=0.\end{cases}
$$

（$\beta_1=0$ 时 $\cos\theta_1=0$，取 $\pi/2$；此时上式给出 $u=CR_1\cos(\sqrt\lambda x)$，与「$u'(0)=0$ 迫使 $B=0$、$A$ 任意」完全一致。）

**第二步：代入 $x=L$ 的边值条件。** 由 $u(L)=CR_1\sin(\sqrt\lambda L+\theta_1)$、$u'(L)=CR_1\sqrt\lambda\cos(\sqrt\lambda L+\theta_1)$，

$$
\alpha_2\sqrt\lambda\cos\bigl(\sqrt\lambda L+\theta_1\bigr)+\beta_2\sin\bigl(\sqrt\lambda L+\theta_1\bigr)=0 .
$$

**两处自变量都是 $\sqrt\lambda L+\theta_1$**：$u$ 与 $u'$ 只差一个 $\pi/2$ 的相位，$x=L$ 处的相位是同一个。合并两个三角项：

$$
R_2\sin\bigl(\sqrt\lambda L+\theta_1+\theta_2\bigr)=0,\qquad R_2=\sqrt{\alpha_2^2\lambda+\beta_2^2},
$$

其中

$$
\sin\theta_2=\frac{\alpha_2\sqrt\lambda}{R_2},\qquad\cos\theta_2=\frac{\beta_2}{R_2},
\qquad
\theta_2=\begin{cases}\arctan\left(\dfrac{\alpha_2\sqrt\lambda}{\beta_2}\right), & \beta_2>0,\\[8pt] \dfrac{\pi}{2}, & \beta_2=0.\end{cases}
$$

**第三步：解相位方程。**

$$
\sqrt\lambda\,L+\theta_1(\lambda)+\theta_2(\lambda)=n\pi,\qquad n=0,1,2,\dots\tag{15.2}
$$

$\theta_i\in[0,\pi/2]$，故 $\sqrt\lambda L\in[n\pi-\pi,\,n\pi]$。

**三类边值的显式谱**

1. **Dirichlet**（$\alpha_1=\alpha_2=0$）：$\theta_1=\theta_2=0$，(15.2) 给出 $\sqrt\lambda=\dfrac{n\pi}{L}$，

$$
\lambda_n=\left(\frac{n\pi}{L}\right)^2,\qquad n=1,2,3,\dots
$$

2. **Neumann**（$\beta_1=\beta_2=0$）：$\theta_1=\theta_2=\dfrac\pi2$，(15.2) 给出 $\sqrt\lambda\,L=(n-1)\pi$，

$$
\lambda_n=\left(\frac{(n-1)\pi}{L}\right)^2,\qquad n=1,2,3,\dots
$$

3. **Robin**（$\beta_i>0$）：

$$
\sqrt\lambda\,L+\arctan\left(\frac{\alpha_1\sqrt\lambda}{\beta_1}\right)+\arctan\left(\frac{\alpha_2\sqrt\lambda}{\beta_2}\right)=n\pi .
$$

左端是 $\sqrt\lambda$ 的严格增函数（三项都增），故每个 $n$ 至多一个解；又 $\sqrt\lambda\to\infty$ 时左端 $\to\infty$，故解存在。于是谱离散且趋于 $+\infty$。

> **要点（$n$ 从几开始 —— 讲稿没交代，这是最容易错的一处）**
> 相位方程 (15.2) 在 $n=0$ 时的形式解是 $\sqrt\lambda=0$（此时 $\theta_1=\theta_2=0$ 只有 Dirichlet 成立）。但 $\lambda=0$ 在 Dirichlet 下只对应零解（$u\equiv$ 常数且 $u(0)=0$），**不是特征值**，必须丢掉。所以 Dirichlet 从 $n=1$ 起、$\lambda_1=(\pi/L)^2>0$；Robin（$\sigma>0$）同样从 $n=1$ 起（$\lambda=0$ 时 $u\equiv c$ 与 $\sigma u|_{\partial}=0$ 迫使 $c=0$）；只有 Neumann 的 $n=1$（对应 $\sqrt\lambda\,L=0$）合法，给出 $\lambda_1=0$ 与常函数。这与定理 15.2.1 的「$\lambda_1=0\iff\beta_1=\beta_2=0$」完全吻合 —— 若忘了丢 $n=0$，Dirichlet 会凭空多出一个 $\lambda=0$。

---

## 15.4 展开定理与 Fourier 级数

**定理 15.4.1**　取 §15.2 的 $\{u_k\}$ 并**归一化**为 $\lVert u_k\rVert_{L^2}=1$、$\int_0^Lu_iu_j=0$（$i\ne j$）。则

1. 对任意 $f\in L^2((0,L))$，

$$
f(x)=\sum_{j=1}^{\infty}a_ju_j(x),\qquad
\left\lVert f-\sum_{j=1}^Na_ju_j\right\rVert_{L^2}\xrightarrow[N\to\infty]{}0 ,
$$

其中 $a_j=\int_0^Lfu_j$。

2. 若 $f\in C^2([0,L])$ 且满足 (15.1) 的边值条件，则 $\sum_j a_ju_j$ 绝对收敛且一致收敛到 $f$。

**补充（具体写出来就是熟悉的 Fourier 级数）**　把三类边值的特征函数显式写出（已归一化）：

| 边值 | 特征函数 | 谱 | 展开名 |
|---|---|---|---|
| Dirichlet | $u_k=\sqrt{2/L}\,\sin\dfrac{k\pi x}{L}$，$k\ge1$ | $(k\pi/L)^2$ | 正弦级数 |
| Neumann | $u_0=1/\sqrt L$；$u_k=\sqrt{2/L}\,\cos\dfrac{k\pi x}{L}$，$k\ge1$ | $0$；$(k\pi/L)^2$ | 余弦级数 |
| Robin | 由相位方程定出的三角函数（相位 $\theta_1$） | 超越方程 (15.2) 的解 | —— |

于是定理 15.4.1 就是「$L^2$ 函数的正弦/余弦展开」。**注意 Neumann 多出一个常数项** $a_0=L^{-1/2}\int f$，它正是 $\lambda=0$ 那个特征方向；这也是为什么 Neumann 的级数在端点处没有「$f(0)=0$」的强制要求，而 Dirichlet 的级数即使 $f$ 连续、只要 $f(0)\ne0$，在端点附近也不会收敛到 $f$ —— 定理 15.4.1 第 2 条要求 $f$ **满足边值条件**，正是为了排除这种情形。

> **要点（第 2 条里「满足边值条件」不能省）**　取 $L=1$、$f\equiv1$、Dirichlet 边值。它的正弦级数在 $x=0$ 处每项都为 $0$，级数只能收敛到 $0\ne f(0)=1$；在 $(0,1)$ 内部也只在 $L^2$ 意义下收敛，端点附近更不可能一致收敛到 $f$。所以相容性要求是实质的，不是技术性的。

---

## 订正

1. **代入右端点边值条件时，一处相位符号写错。** 讲稿写 $\alpha_2\sqrt\lambda\cos(\sqrt\lambda L+\theta_1)+\beta_2\sin(\sqrt\lambda L+\theta_2)=0$，第二项应为 $\sin(\sqrt\lambda L+\theta_1)$。理由：$u(L)$ 与 $u'(L)$ 都由**同一个** $\sin(\sqrt\lambda L+\theta_1)$ 及其导数给出，$\theta_2$ 是**下一步**（合并两个三角项）才出现的相位。照原式写，合并出的相位就不是 $\sqrt\lambda L+\theta_1+\theta_2$，与紧随其后的结论矛盾。
2. **第七周把一致收敛的条件简写成 $f\in C^2([0,L])$，漏了「且 $f$ 满足 (15.1) 的边值条件」。** 第六周 §3.3 第 3 条写全了。反例见 §15.4 的要点：$f\equiv1$ 在 Dirichlet 边值下就不满足相容性。这个遗漏会让人以为「只要 $f$ 光滑就一致收敛」，进而误判端点行为。
3. **Neumann 的谱 $\lambda_n=((n-1)\pi/L)^2$ 没注明 $n\ge1$。** 若把 $n$ 从 $0$ 起算，会得到 $(-\pi/L)^2$ 与 $n=1$ 重复，看似无碍，实则把「按重数排列」的编号打乱：$\lambda_1=0$ 应是单重（常数），照 $n=0,1$ 数会算成双重。正确写法 $n=1,2,3,\dots$。
4. **定理 15.2.1 第 2 条说 $\{u_k\}$ 构成「完备标准正交基」，但 §15.3 算出的特征函数没有归一化。** 归一化常数必须显式写出（见 §15.4 的表），否则展开系数 $a_k=\int fu_k$ 是错的（应除以 $\lVert u_k\rVert^2$）。第七周那版补了 $\lVert u_k\rVert=1$ 与 $\int u_iu_j=0$，第六周那版没有。
5. **两讲重复陈述同一个定理，且措辞、编号含义不一致。** 第六周版本把「$\lambda_1=0$ 的判别」放在「特别地」里，第七周把它并进主句；两处对 $n$ 的起算也未统一（第六周未提，第七周用 $n$ 通吃 $0,1,2,\dots$）。本篇统一为 $n\ge1$ 并在 §15.3 的要点里说明 $n=0$ 分支为何要丢。
6. **§15.3 里 $u=C[\alpha_1\sqrt\lambda\cos(\sqrt\lambda x)+\beta_1\sin(\sqrt\lambda x)]$ 这一步是「除以 $\beta_1$」得到的，$\beta_1=0$ 时该推导失效。** 讲稿没交代。正确做法：$\beta_1=0$ 时由 $-\alpha_1\sqrt\lambda\,B=0$ 直接得 $B=0$、$u=A\cos(\sqrt\lambda x)$，与通式取 $\beta_1=0$ 的结果 $C\alpha_1\sqrt\lambda\cos(\sqrt\lambda x)$ 一致。所以通式在 $\beta_1=0$ 时**仍然成立**，只是不能靠「除」得到，要靠直接代入。

---

系列：[[14-谱理论·鼓面的数学建模与Laplace特征值]] · [[16-特征值的变分刻画与Poincaré不等式]] · [[12-习题课]]
