# 10 球与半空间上的 Green 函数

> 系列：[[09-Green函数的性质]] · [[11-Poisson积分公式与Dirichlet问题的可解性]] · [[12-习题课]]

[[08-Green恒等式与Green表示公式]] 把解 Dirichlet 问题归结为：造出

$$
G(x,y)=\Gamma(x-y)-\gamma(x,y),\qquad
\begin{cases}
\Delta_y\gamma(x,y)=0, & y\in\Omega,\\[4pt]
\gamma(x,y)=\Gamma(x-y), & y\in\partial\Omega .
\end{cases}
$$

一般区域上 $\gamma$ 写不出来；但**球**和**半空间**上因为对称性足够高，可以用「对称点」把 $\gamma$ 显式地写出来。本篇就做这件事。

---

## 10.1 球 $B_R$ 上的 Green 函数

**定理**　设 $B_R(0)\subseteq\mathbb R^n$，则对于 $\forall x\in B_R\setminus\{0\}$、$y\in B_R$、$y\ne x$，

$$
G(x,y)=\begin{cases}
\dfrac{1}{2\pi}\left(\log\dfrac{1}{|x-y|}+\log\left|\dfrac{R}{|x|}x-\dfrac{|x|}{R}y\right|\right), & n=2,\\[14pt]
\dfrac{1}{(n-2)|S^{n-1}|}\left(\dfrac{1}{|x-y|^{n-2}}-\dfrac{1}{\left|\dfrac{R}{|x|}x-\dfrac{|x|}{R}y\right|^{n-2}}\right), & n\ge3 .
\end{cases}
$$

**证明**　关键是找到 $\gamma(x,y)$，使得

$$
\begin{cases}
\Delta_y\gamma(x,y)=0, & y\in B_R,\\[4pt]
\gamma(x,y)=\Gamma(x-y), & y\in\partial B_R .
\end{cases}
$$

**$\Gamma(x-y)$ 与我们想要的 $\gamma(x,y)$ 差在哪？** 差在 $\Gamma(x-y)$ 有奇点 $y=x$，而 $\gamma$ 要求在整个 $B_R$（含球心附近）调和。

**利用球面对称：对称点。** 任取 $x\in B_R\setminus\{0\}$，令 $x^*$ 为 $x$ 关于球面 $\partial B_R$ 的**对称点**：

$$
|x^*|\cdot|x|=R^2,\qquad\text{即}\qquad
x^*=\left(\frac{R^2}{|x|}\right)\frac{x}{|x|}=\frac{R^2}{|x|^2}\,x .
$$

（注意 $|x^*|>R$，故 $x^*\notin B_R$——这正是「$\Gamma(x^*-y)$ 在 $B_R$ 内没有奇点」的原因。）

**相似三角形给出边界上的距离关系。** 由于

$$
\frac{|x^*|}{R}=\frac{R}{|x|}
$$

且 $\triangle Oxy$ 与 $\triangle Oyx^*$ 有公共角 $\angle xOy$、夹这个角的两边成同一比例 $|Ox|/|Oy|=|Oy|/|Ox^*|=|x|/R$，故二者相似。因此在对应边上

$$
\frac{|x^*-y|}{|x-y|}=\frac{|x^*|}{R}=\frac{R}{|x|},
\qquad\text{即}\qquad
|x^*-y|=\frac{R}{|x|}|x-y|,
$$

对任一点 $y\in\partial B_R$ 成立。

**分两种情形构造 $\gamma$。**

**$n\ge3$：用齐次性。** 此时 $\Gamma$ 满足

$$
\Gamma(\lambda r)=\lambda^{-(n-2)}\Gamma(r),
$$

于是当 $y\in\partial B_R$ 时

$$
\Gamma(x^*-y)=\left(\frac{|x|}{R}\right)^{n-2}\Gamma(x-y).
$$

两边同乘 $(R/|x|)^{n-2}$ 即得边界条件。而 $x^*\notin B_R$，故

$$
\gamma(x,y):=\left(\frac{R}{|x|}\right)^{n-2}\Gamma(x^*-y),
$$

**关于 $y$ 在 $B_R$ 内调和**（因为 $\Gamma$ 在 $\mathbb R^n\setminus\{0\}$ 上调和、$x^*-y\ne0$ 于 $B_R$），且在 $\partial B_R$ 上等于 $\Gamma(x-y)$。

**$n=2$：用 $\log$ 的加法形式。** 此时

$$
\Gamma(\lambda r)=-\frac{1}{2\pi}\log\lambda+\Gamma(r),
$$

于是当 $y\in\partial B_R$ 时 $\Gamma(x^*-y)=-\dfrac{1}{2\pi}\log\dfrac{R}{|x|}+\Gamma(x-y)$。故令

$$
\gamma(x,y):=\Gamma(x^*-y)+\frac{1}{2\pi}\log\frac{R}{|x|},
$$

由 $x^*\notin B_R$ 与边界上的距离关系知 $\gamma(x,\cdot)$ 在 $B_R$ 内调和，且在 $\partial B_R$ 上等于 $\Gamma(x-y)$。

**化简。** 代入 $x^*=\dfrac{R^2}{|x|^2}x$。注意

$$
\frac{R}{|x|}x-\frac{|x|}{R}y=\frac{|x|}{R}\left(x^*-y\right),
$$

故 $\left|\dfrac{R}{|x|}x-\dfrac{|x|}{R}y\right|=\dfrac{|x|}{R}|x^*-y|$。于是

$$
\gamma(x,y)=\Gamma\!\left(\frac{|x|}{R}\,|x^*-y|\right)
=\frac{1}{(n-2)|S^{n-1}|}\cdot\frac{1}{\left|\dfrac{R}{|x|}x-\dfrac{|x|}{R}y\right|^{n-2}}\quad(n\ge3),
$$

> 对 $n=2$ 同理：$\Gamma\!\left(\frac{|x|}{R}|x^*-y|\right)=-\frac{1}{2\pi}\log\frac{|x|}{R}-\frac{1}{2\pi}\log|x^*-y|=\Gamma(x^*-y)+\frac{1}{2\pi}\log\frac{R}{|x|}$，与上面一致。

把它代回 $G=\Gamma(x-y)-\gamma(x,y)$，即得定理中的两式。$\square$

> **要点**　两式可以**统一写成**

$$
G(x,y)=\Gamma(x-y)-\Gamma(z),\qquad z=\frac{R}{|x|}x-\frac{|x|}{R}y,
$$

> 也就是说：$n\ge2$ 一律有 $G=\Gamma(x-y)-\Gamma(\frac{|x|}{R}|x^*-y|)$，只是 $n=2$ 时 $\Gamma$ 的齐次性是「加法形式」而非「乘法形式」。这一步「$n=2$ 与 $n\ge3$ 的差别只在 $\Gamma$ 的齐次性」是本篇最值得记住的结构。

### 10.1.1 $x=0$ 处的延拓与对称性

**$x=0$ 的延拓**　当 $x\to0$ 时 $|z|\to1$、$|x-y|\to|y|$，故上式在 $x=0$ 处延拓为

$$
G(0,y)=\Gamma(y)-\Gamma(1)=\begin{cases}
-\dfrac{1}{2\pi}\log|y|, & n=2,\\[8pt]
\dfrac{|y|^{2-n}-1}{(n-2)|S^{n-1}|}, & n\ge3,
\end{cases}
\qquad0<|y|\le R
$$

> （取 $R=1$ 写出；一般 $R$ 时把 $\Gamma(1)$ 换成 $\Gamma(R)$ 即可。）

**对称性**　注意

$$
\left|\frac{R}{|x|}x-\frac{|x|}{R}y\right|^2=R^2-2\langle x,y\rangle+\frac{|x|^2|y|^2}{R^2}
$$

（对 $x$ 与 $y$ 完全对称），故 $G(x,y)=G(y,x)$——与 [[09-Green函数的性质]] 性质 3 一致。

> **补充（$R=1$ 时正性的一个一步证明）**　对 $R=1$，把上式与 $|x-y|^2=|x|^2-2\langle x,y\rangle+|y|^2$ 相减：

$$
\left|\frac{x}{|x|}-|x|y\right|^2-|x-y|^2=\bigl(1-|x|^2\bigr)\bigl(1-|y|^2\bigr)\ge0\qquad(x,y\in B_1).
$$

> 因为 $\Gamma$ 关于半径**单调递减**，故 $\Gamma(|z|)\le\Gamma(x-y)$，即

$$
G(x,y)=\Gamma(x-y)-\Gamma(|z|)\ge0,\qquad x,y\in B_1 .
$$

> 这给了 [[09-Green函数的性质]] 性质 2 的正性在球上的一个独立而简洁的证明（不依赖极值原理）。这条恒等式在 [[11-Poisson积分公式与Dirichlet问题的可解性]] 处理 $-\Delta u=f$、$u|_{\partial B_1}=0$ 时还要用到。

---

## 10.2 半空间 $\mathbb R^n_+$ 上的 Green 函数

用**镜像法**。对于 $x\in\mathbb R^n_+$，仍要求 $\gamma(x,y)$ 满足

$$
\begin{cases}
\Delta_y\gamma(x,y)=0, & y\in\mathbb R^n_+,\\[4pt]
\gamma(x,y)=\Gamma(x-y), & y\in\partial\mathbb R^n_+=\{y_n=0\} .
\end{cases}
$$

取 $x^*$ 为 $x$ 关于超平面 $\mathbb R^{n-1}$ 的对称点：

$$
x^*=(x_1,\dots,x_{n-1},-x_n),
$$

则当 $y_n=0$ 时（此时 $y=x^*$ 的镜像也有 $y_n=0$）

$$
|x^*-y|=|x-y| .
$$

又 $x^*\notin\mathbb R^n_+$，故 $y\mapsto\Gamma(x^*-y)$ 在整个 $\mathbb R^n_+$ 上调和。因此直接取

$$
\gamma(x,y)=\Gamma(x^*-y)
$$

即可。于是

$$
G(x,y)=\Gamma(x-y)-\Gamma(x^*-y),\qquad x^*=(x_1,\dots,x_{n-1},-x_n).
$$

> **要点**　半空间与球的做法是**同一件事的两个化身**：
>
> | 区域 | 「反射」 | 保距的对称 |
> |---|---|---|
> | 半空间 $\mathbb R^n_+$ | 关于超平面反射 $x\mapsto(x_1,\dots,x_{n-1},-x_n)$ | 反射是等距，故 $\lvert x^*-y\rvert=\lvert x-y\rvert$ 于边界 |
> | 球 $B_R$ | 关于球面反演 $x\mapsto\dfrac{R^2}{\lvert x\rvert^2}x$ | 反演不是等距，只在**球面上**满足 $\lvert x^*-y\rvert=\dfrac{R}{\lvert x\rvert}\lvert x-y\rvert$ |
>
> 「边界上距离成定比」是能写出 $\gamma$ 的全部原因。对一般区域，边界上这个比例不再是常数，所以 $\gamma$ 没有显式解——只能靠 [[11-Poisson积分公式与Dirichlet问题的可解性]] 里那种「先假定有解再验证」的办法论证存在性。

> **补充（相似三角形的严格论证）**　讲稿直接断言 $\triangle Oxy$ 与 $\triangle Oyx^*$ 相似。严格说：两三角形有公共角 $\angle xOy$，且夹该角的两边

$$
\frac{|Ox|}{|Oy|}=\frac{|x|}{R},\qquad
\frac{|Oy|}{|Ox^*|}=\frac{R}{R^2/|x|}=\frac{|x|}{R}
$$

> 成同一比例，故由「两边成比例且夹角相等」得相似，相似比 $|x|/R$。对应边 $xy\leftrightarrow yx^*$，故

$$
\frac{|x^*-y|}{|x-y|}=\frac{|Ox^*|}{|Oy|}=\frac{R^2/|x|}{R}=\frac{R}{|x|}.
$$

> 注意这里必须用 $|x^*|=R^2/|x|$（即 $|x^*||x|=R^2$），这是对称点的定义。讲稿把这个等式写成 $|x^*|/R=R/|x|$，是同一件事。

---

## 订正

1. **讲稿 §5 的编号「（1）定理」悬空。** 该节标题是「在 $B_R$ 和 $\mathbb R^n_+$ 上求解 Green 函数」，开头标了「（1）定理」，但后面讲半空间时并未出现「（2）」，而是直接写「类似地，我们可得到 $\mathbb R^n_+$ 上的 Green 函数」。本篇把两者分别列为 §10.1 与 §10.2。
2. **$\Gamma$ 的齐次性，两种形式不能混用。** 讲稿两处分别写出
   $\Gamma(\lambda r)=\lambda^{-(n-2)}\Gamma(r)$（$n\ge3$）与 $\Gamma(\lambda r)=-\frac{1}{2\pi}\log\lambda+\Gamma(r)$（$n=2$）。第二个式子常被误读成「$n=2$ 时也有幂次形式」，从而在 $n=2$ 的 $\gamma$ 里漏掉 $+\frac{1}{2\pi}\log\frac{R}{|x|}$ 这一项（漏掉后 $\gamma$ 在球面上就不等于 $\Gamma(x-y)$，$G$ 的零边值不成立）。本篇已在证明中把两种情形分开写。
3. **相似三角形一句未给理由。** 讲稿写「由 $|x^*|/R=R/|x|$ 和 $\triangle Oxy$ 与 $\triangle Oyx^*$ 为相似三角形」，但相似性本身需要「公共角 + 夹边成比例」才成立。本篇已在**补充**中补出。
4. **$x=0$ 的延拓只给了 $R=1$ 的具体值。** 讲稿（第四周）在球上写「当 $x=0$ 时，上式延拓为 $G(0,y)=\Gamma(y)-\Gamma(1)$」，用的是 $R=1$；第五周又用 $G(0,y)=\Gamma(y)-\Gamma(1)$ 验证 Poisson 核的边界行为，也是 $R=1$。一般 $R$ 时的对应式是 $G(0,y)=\Gamma(y)-\Gamma(R)$，本篇已注明。

---

系列：[[09-Green函数的性质]] · [[11-Poisson积分公式与Dirichlet问题的可解性]] · [[12-习题课]]
