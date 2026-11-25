# 11 Poisson 积分公式与 Dirichlet 问题的可解性

> 系列：[[10-球与半空间上的Green函数]] · [[09-Green函数的性质]] · [[12-习题课]]

[[10-球与半空间上的Green函数]] 造出了球上的 $G$。本篇做两件事：

1. 算出 $\dfrac{\partial G}{\partial\nu_y}$ 在 $\partial B_1$ 上的显式式（**Poisson 核**）；
2. 用它**验证**球上 Dirichlet 问题的可解性——即证明 [[08-Green恒等式与Green表示公式]] 中那个「形式上的求解公式」确实给出解。

讲义里这两件事分散在两讲（先算核并验证 $-\Delta u=0$、$u|_{\partial B_1}=g$；再补上 $-\Delta u=f$、$u|_{\partial B_1}=0$ 的情形），本篇按知识合为一处。

---

## 11.1 球上的 Poisson 核

在 $B_1$ 上（即 $R=1$），由 [[10-球与半空间上的Green函数]] §10.1，

$$
G(x,y)=\Gamma(x-y)-\Gamma(z),\qquad
z=\frac{x}{|x|}-|x|y,\qquad
|z|^2=1-2\langle x,y\rangle+|x|^2|y|^2 .
$$

**命题**　对 $x\in B_1$、$y\in\partial B_1$，

$$
\frac{\partial G}{\partial\nu_y}(x,y)=-\frac{1}{|S^{n-1}|}\cdot\frac{1-|x|^2}{|x-y|^n}<0 .
$$

**证明**　任取 $y\in\partial B_1$，$\nu_y=\dfrac{y}{|y|}=y$。

**$n\ge3$ 的情形。** 把 $G$ 写成 $\Gamma(x-y)-\left(\dfrac{1}{|x|}\right)^{n-2}\Gamma(x^*-y)$（$x^*=x/|x|^2$），用

$$
\nabla_y\Gamma(x-y)=\frac{x-y}{|S^{n-1}||x-y|^n},\qquad
\nabla_y\Gamma(x^*-y)=\frac{x^*-y}{|S^{n-1}||x^*-y|^n},
$$

得

$$
\frac{\partial G}{\partial\nu_y}(x,y)
=\frac{1}{|S^{n-1}|}\left[-\frac{y-x}{|x-y|^n}+\frac{y-x^*}{|x|^{n-2}|x^*-y|^n}\right]\cdot\frac{y}{|y|},
$$

其中第一项已把 $(x-y)$ 写成 $-(y-x)$。

**两种情形都可以化到同一个式子。** 在 $\partial B_1$ 上由相似三角形关系 $|x^*-y|=\dfrac{|x-y|}{|x|}$（见 [[10-球与半空间上的Green函数]]）得

$$
|x|^{n-2}|x^*-y|^n=|x|^{n-2}\cdot\frac{|x-y|^n}{|x|^n}=\frac{|x-y|^n}{|x|^2},
$$

故

$$
\frac{y-x^*}{|x|^{n-2}|x^*-y|^n}
=\frac{|x|^2(y-x^*)}{|x-y|^n}
=\frac{|x|^2y-x}{|x-y|^n}\qquad(\text{用了 }|x|^2x^*=x),
$$

于是

$$
\frac{\partial G}{\partial\nu_y}(x,y)
=\frac{1}{|S^{n-1}||x-y|^n}\bigl[-(y-x)+(|x|^2y-x)\bigr]\cdot y
=\frac{(|x|^2-1)\langle y,y\rangle}{|S^{n-1}||x-y|^n}
=-\frac{1-|x|^2}{|S^{n-1}||x-y|^n}.
$$

（$n=2$ 时把 $|S^{n-1}|$ 换成 $2\pi$，同一计算给出同一个式子。）又因 $|y|=1$ 使 $\langle y,y\rangle=1$。$\square$

**$x=0$ 的相容性**　由 $G(0,y)=\Gamma(y)-\Gamma(1)$（[[10-球与半空间上的Green函数]] §10.1.1）直接算：$\nabla_yG(0,y)=\nabla\Gamma(y)=-\dfrac{y}{|S^{n-1}|}$（$|y|=1$），故

$$
\frac{\partial G}{\partial\nu_y}(0,y)=\left\langle-\frac{y}{|S^{n-1}|},y\right\rangle=-\frac{1}{|S^{n-1}|},
$$

与上面的公式在 $x=0$ 处的值 $-\dfrac{1-0}{|S^{n-1}|\cdot1}=-\dfrac{1}{|S^{n-1}|}$ 一致。**所以该公式在 $x=0$ 也成立。**

由 [[09-Green函数的性质]] 性质 4，

$$
\int_{\partial B_1}\frac{\partial G}{\partial\nu_y}(x,y)\,d\sigma_y=-1,\qquad\forall x\in B_1 .
$$

> **要点**　$\dfrac{\partial G}{\partial\nu_y}$ **严格为负**（因为 $|x|<1$）且**总积分为 $-1$**。两条合起来说明 $-\dfrac{\partial G}{\partial\nu_y}\ge0$ 且积分为 $1$——即 $-\dfrac{\partial G}{\partial\nu_y}$ 是 $\partial B_1$ 上的一个**概率核**。这是「边界值的加权平均」这一直觉的来源，也是下面验证 $u|_{\partial B_1}=g$ 的关键。

**Poisson 积分公式**　把核代回 (8.1)，得球上调和函数的显式表示

$$
u(x)=\frac{1}{|S^{n-1}|}\int_{\partial B_1}\frac{1-|x|^2}{|x-y|^n}\,u(y)\,d\sigma_y\qquad(x\in B_1,\ \Delta u=0).
$$

---

## 11.2 定理：$\Delta u=0$、$u|_{\partial B_1}=g$ 可解

**定理**　令

$$
u(x)=-\int_{\partial B_1}g(y)\frac{\partial G}{\partial\nu_y}(x,y)\,d\sigma_y,\qquad x\in B_1(0),
$$

其中 $G(x,y)$ 为球 $B_1$ 的 Green 函数，$g\in C(\partial B_1)$。则 $u$ 是下述调和方程 Dirichlet 问题的解：

$$
\begin{cases}
\Delta u=0, & x\in B_1,\\
u=g, & x\in\partial B_1 .
\end{cases}
$$

**证明**　分两步。

**第一步：$u\in C^2(B_1)$ 且 $\Delta u=0$。** 设

$$
u(x)=-\int_{\partial B_1}g(y)\frac{\partial G}{\partial\nu_y}(x,y)\,d\sigma_y .
$$

当 $x=0$ 时取 $G(0,y)=\Gamma(y)-\Gamma(1)$（$y\ne0$），即前式在 $x=0$ 处的光滑延拓。对任意 $K\Subset B_1$，$G$ 在 $K$ 与 $\partial B_1$ 的一个邻域的乘积上光滑，且关于 $x$ 调和，故核及其关于 $x$ 的各阶导数在 $K\times\partial B_1$ 上连续有界，**可以在积分号下求导**：

$$
\Delta_xu(x)=-\int_{\partial B_1}g(y)\frac{\partial}{\partial\nu_y}\bigl(\Delta_xG(x,y)\bigr)d\sigma_y=0,
$$

且 $u\in C^2(B_1)$。（$\Delta_xG=0$ 于 $x\ne y$，且在 $K\times\partial B_1$ 上 $x\ne y$。）

**第二步：验证边界值 $u(x)\to g(x_0)$（$x\to x_0$，$x_0\in\partial B_1$）。** 由 §11.1 的核公式与通量等式（$\int_{\partial B_1}\partial G/\partial\nu_y\,d\sigma_y=-1$），

$$
u(x)-g(x_0)=u(x)+g(x_0)\int_{\partial B_1}\frac{\partial G}{\partial\nu_y}(x,y)\,d\sigma_y
=-\int_{\partial B_1}\bigl[g(y)-g(x_0)\bigr]\frac{\partial G}{\partial\nu_y}(x,y)\,d\sigma_y .
$$

**目标：证明 $x\to x_0$ 时 $\bigl|u(x)-g(x_0)\bigr|\to0$。** 把积分按 $|y-x_0|$ 是否超过 $|x-x_0|^{1/2}$ 拆为**近区与远区**：

$$
u(x)-g(x_0)=-\int_{\partial B_1\cap\{|y-x_0|\le|x-x_0|^{1/2}\}}\bigl[g(y)-g(x_0)\bigr]\frac{\partial G}{\partial\nu_y}\,d\sigma_y
-\int_{\partial B_1\cap\{|y-x_0|>|x-x_0|^{1/2}\}}\bigl[g(y)-g(x_0)\bigr]\frac{\partial G}{\partial\nu_y}\,d\sigma_y
=:I+II .
$$

**近区 $I$。** 由 $\dfrac{\partial G}{\partial\nu_y}<0$ 与通量等式 $-1$ 及 $g$ 的连续性，

$$
|I(x)|\le\sup_{|y-x_0|\le|x-x_0|^{1/2}}\bigl|g(y)-g(x_0)\bigr|\cdot\left(-\int_{\partial B_1}\frac{\partial G}{\partial\nu_y}\,d\sigma_y\right)
=\sup_{|y-x_0|\le|x-x_0|^{1/2}}\bigl|g(y)-g(x_0)\bigr|\longrightarrow0 .
$$

**远区 $II$。** 当 $0<|x-x_0|<\dfrac14$ 时，在远区上有

$$
|x-y|\ge\bigl||y-x_0|-|x-x_0|\bigr|\ge\frac12|y-x_0|>\frac12|x-x_0|^{1/2},
$$

因而

$$
-\frac{\partial G}{\partial\nu_y}(x,y)=\frac{1}{|S^{n-1}|}\frac{1-|x|^2}{|x-y|^n}
\le\frac{2^n}{|S^{n-1}|}\cdot\frac{1-|x|^2}{|y-x_0|^n}
\le\frac{2^{n+1}}{|S^{n-1}|}\cdot\frac{|x-x_0|}{|y-x_0|^n},
$$

最后一步用 $1-|x|^2\le2\bigl|(1-|x|)\bigr|\le2|x-x_0|$。于是

$$
|II(x)|\le C_n\|g\|_{L^\infty(\partial B_1)}\,|x-x_0|\int_{\partial B_1\cap\{|y-x_0|>|x-x_0|^{1/2}\}}\frac{d\sigma_y}{|y-x_0|^n}.
$$

**估最后这个积分。** 取 $r=|x-x_0|^{1/2}$。用「球面帽的测度控制」：存在 $C_n$ 使对一切 $0<s\le2$，

$$
\bigl|\partial B_1\cap B_s(x_0)\bigr|\le C_ns^{\,n-1}
$$

（$s\le1$ 时，把球面帽投影到 $x_0$ 处的切平面，球面面积元至多为投影面积元的两倍；$s\ge1$ 时用整球面的面积即可）。再按二进制环带分解：

$$
\int_{\partial B_1\cap\{|y-x_0|>r\}}\frac{d\sigma_y}{|y-x_0|^n}
\le\sum_{j=0}^\infty\frac{\bigl|\partial B_1\cap B_{2^{j+1}r}(x_0)\bigr|}{(2^jr)^n}
\le\frac{C_n}{r}\sum_{j=0}^\infty2^{-j}\le\frac{C_n}{r}=\frac{C_n}{|x-x_0|^{1/2}} .
$$

代回得

$$
|II(x)|\le C_n\|g\|_{L^\infty(\partial B_1)}\,|x-x_0|^{1/2}\longrightarrow0\qquad(x\to x_0).
$$

因此 $\lim_{x\to x_0}u(x)=g(x_0)$。令 $u|_{\partial B_1}=g$，则 $u\in C^2(B_1)\cap C(\overline{B_1})$，满足所述的 Dirichlet 问题。$\square$

> **补充（另一种拆法）**　讲义在另一处把远区拆成 $|y-x_0|>|x-x_0|^{1/(2n)}$，此时

$$
|x-y|\ge|x-x_0|^{1/(2n)}-|x-x_0|\ge\frac14|x-x_0|^{1/(2n)},
$$

> 于是 $-\partial G/\partial\nu_y\le C_n|x-x_0|^{1/2}$ 是**逐点**的界，远区的测度又是常数，直接得 $|II|\le C_n\|g\|_\infty|x-x_0|^{1/2}$。这种拆法更省事（不需要二进制环带），代价是 $1/(2n)$ 的指数看着突兀。两种拆法都对，本篇给出的是「近区用连续性、远区用核的衰减」这一条主线。

> **要点**　第一步与第二步的技术含量不同：**$\Delta u=0$ 靠「在积分号下求导」（核在 $K\times\partial B_1$ 上光滑）**，**边界值靠「核的正则化性质」（非负、积分为 $1$、质量集中到边界点）**。后者与 [[07-基本解与广义函数]] 里 $q_\varepsilon$ 的「近似恒等核」完全相同——只是从球面上的核换成了边界上的核。

---

## 11.3 情形 $\Delta u=f$、$u|_{\partial B_1}=0$

**定理**　设 $f\in C^1(B_1)$，令

$$
u(x)=\int_{B_1}f(y)G(x,y)\,dy,\qquad x\in B_1 .
$$

则 $u\in C^2(B_1)\cap C(\overline{B_1})$，且

$$
\begin{cases}
-\Delta u=f, & x\in B_1,\\
u=0, & x\in\partial B_1 .
\end{cases}
$$

**证明**　分三步。

**第一步：$G\ge0$ 于 $B_1\times B_1$。** 由 [[10-球与半空间上的Green函数]] §10.1.1 的恒等式

$$
1-2\langle x,y\rangle+|x|^2|y|^2-|x-y|^2=\bigl(1-|x|^2\bigr)\bigl(1-|y|^2\bigr)\ge0
$$

及 $\Gamma$ 关于半径的单调递减性知 $G(x,y)=\Gamma(x-y)-\Gamma(|z|)\ge0$。

**第二步：边界上 $u\to0$。** 令

$$
w(x)=\frac{1-|x|^2}{2n},
$$

则 $\Delta w=-\dfrac{2n}{2n}=-1$，即 $-\Delta w=1$，且 $w|_{\partial B_1}=0$。**把 Green 表示公式 (8.1) 用于 $w$**（此时 $-\Delta w=1$、$w|_{\partial B_1}=0$，含 $\partial w/\partial\nu$ 的项被 $G|_{\partial B_1}=0$ 消掉）：

$$
w(x)=\int_{B_1}1\cdot G(x,y)\,dy .
$$

于是

$$
\int_{B_1}G(x,y)\,dy=\frac{1-|x|^2}{2n},
$$

从而

$$
|u(x)|\le\|f\|_{L^\infty(B_1)}\int_{B_1}G(x,y)\,dy\le\|f\|_{L^\infty(B_1)}\frac{1-|x|^2}{2n}\longrightarrow0\qquad(x\to\partial B_1).
$$

所以把 $u$ 在边界上定义为 $0$ 后，$u$ 在边界处连续。

**第三步：验证 $-\Delta_xu(x)=f(x)$。** 把 $G=\Gamma-\gamma$ 拆开：

$$
u(x)=\int_{B_1}f(y)\Gamma(x-y)\,dy-\int_{B_1}f(y)\gamma(x,y)\,dy .
$$

**后一项是调和函数。** 注意

$$
\gamma(x,y)=\Gamma\!\left(|y|\,|x-y^*|\right),\qquad y^*=\frac{y}{|y|^2}\notin B_1
$$

（由 [[10-球与半空间上的Green函数]] 中 $|z|$ 的两条对称表达式可得），故 $x\mapsto\gamma(x,y)$ 在 $B_1$ 上调和。由控制收敛定理（$\gamma$ 在 $\overline{B_1}\times\overline{B_1}$ 上光滑、各阶导数有界）可在积分号下求导，故 $\displaystyle\int_{B_1}f(y)\gamma(x,y)dy$ 调和。

**前一项要用正则化。** $\Gamma(x-y)$ 的二阶导数不能用控制收敛（奇性太强），但一阶导数可积。设

$$
\Gamma_\varepsilon(x)=\begin{cases}
-\dfrac{1}{2\pi}\log\sqrt{|x|^2+\varepsilon^2}, & n=2,\\[10pt]
\dfrac{1}{(n-2)|S^{n-1}|\bigl(|x|^2+\varepsilon^2\bigr)^{(n-2)/2}}, & n\ge3,
\end{cases}
\qquad
q_\varepsilon(x):=-\Delta\Gamma_\varepsilon(x)=\frac{n\varepsilon^2}{|S^{n-1}|\bigl(|x|^2+\varepsilon^2\bigr)^{n/2+1}} .
$$

由 $q_\varepsilon$ 的**近似单位性质**（$\int_{\mathbb R^n}q_\varepsilon=1$、$q_\varepsilon\ge0$、质量集中到原点——即 [[07-基本解与广义函数]] 的命题 1、2，也是第二次作业第 2 题），

$$
-\Delta_x\int_{B_1}f(y)\Gamma_\varepsilon(x-y)\,dy=\int_{B_1}f(y)\,q_\varepsilon(x-y)\,dy\longrightarrow f(x)\qquad(\varepsilon\to0^+),
$$

最后一步用了 $f$ 在 $x$ 处连续。于是

$$
-\Delta_xu(x)=f(x)-0=f(x),\qquad\forall x\in B_1 .
$$

$\square$

> **补充（假设可以放宽）**　第一步的 $f\in C^1$ 其实是用不完的：第二步只需 $f\in L^\infty$，第三步的调和部分只需 $f\in L^1$、近似单位部分只需 $f$ 在 $x$ 处连续。所以 $f\in C(B_1)$ 就够了。讲义取 $C^1$ 是为了叙述简洁。

> **补充（$\int_{B_1}G(x,y)dy$ 这条恒等式的另一种用法）**　它说明 $G(x,\cdot)$ 对 $y$ 的总质量是「$-\Delta w=1$ 的解在 $x$ 处的值」，也就是球上 Green 函数对常数源项的响应。这条恒等式在验证 Poisson 核的归一化、以及 [[12-习题课]] 里的若干练习题中反复出现。

---

## 11.4 合并：球上 Dirichlet 问题的可解性

把 §11.2 与 §11.3 线性叠加：设 $f\in C^1(B_1)$、$g\in C(\partial B_1)$，

$$
u(x)=\int_{B_1}f(y)G(x,y)\,dy-\int_{\partial B_1}g(y)\frac{\partial G}{\partial\nu_y}(x,y)\,d\sigma_y,\tag{11.1}
$$

则 $u\in C^2(B_1)\cap C^0(\overline{B_1})$，且

$$
\begin{cases}
-\Delta u=f, & x\in B_1,\\
u=g, & x\in\partial B_1 .
\end{cases}
$$

**证明**　把问题拆成两个：

$$
D_1:\begin{cases}\Delta u=0, & x\in B_1,\\ u=g, & x\in\partial B_1,\end{cases}
\qquad
D_2:\begin{cases}-\Delta u=f, & x\in B_1,\\ u=0, & x\in\partial B_1 .\end{cases}
$$

由 §11.2，$D_1$ 的解由 (11.1) 中的边界积分项给出；由 §11.3，$D_2$ 的解由 (11.1) 中的体积分项给出。两者相加即为 (11.1)，由线性性它满足原问题。$\square$

> **要点**　至此，「**球上 Dirichlet 问题存在唯一解**」这件事被完整证明了：存在性由 (11.1) 给出（构造性），唯一性由 [[03-调和函数·平均值性质与极值原理]] 推论 1 给出。注意解的存在性论证里**没有用到任何函数空间泛函分析**——它是靠 $G$ 的显式表达式一步步算出来的。这也说明为什么「能把 $G$ 写出来」如此重要。

---

## 订正

1. **两讲把同一件事拆成了两半，并使用了不同的拆分指数。** 先讲只验证了 $D_1$（$\Delta u=0$、$u|_{\partial B_1}=g$），用 $|x-x_0|^{1/2}$ 拆分并对远区做二进制环带估计；后讲验证 $D_1+D_2$，改用 $|x-x_0|^{1/(2n)}$ 拆分。两处都对，但同一条证明出现两种指数会让复习时误以为其中一种错了。本篇按知识合并为一节，主用 $1/2$，并在**补充**里给出 $1/(2n)$ 的拆法。
2. **Poisson 核的推导里，讲稿的分子一度印成 $|x_0|^2-|x|^2$。** 在单位球面上 $|x_0|=1$，故它与 $1-|x|^2$ 是同一个量，**并非错误**；但换到 $B_R$（$R\ne1$）就会不自洽。本篇统一写成 $1-|x|^2$（单位球情形）。
3. **$\partial G/\partial\nu_y$ 在 $x=0$ 处成立，讲稿只说「由 $G(0,y)$ 的表达式知」，未给验算。** 本篇已把 $x=0$ 的验算写出。
4. **§11.3 中「后者由控制收敛定理知为调和函数」一句未点明为什么 $\gamma(\cdot,y)$ 在 $x$ 上调和。** 讲稿只把 $\gamma$ 当成 $y$ 的调和函数（那是它的定义）。要把 $\int f(y)\gamma(x,y)dy$ 说成 $x$ 的调和函数，必须用另一条对称表达式 $\gamma(x,y)=\Gamma(|y||x-y^*|)$、$y^*\notin B_1$。本篇已补。（这一步若不补，读者会以为 $\gamma$ 只在 $y$ 变量上调和。）
5. **拆 $D_1/D_2$ 的理由讲稿未说明。** 两者的技术手段完全不同：$D_1$ 靠核在 $K\times\partial B_1$ 上光滑（可积分号下求导），$D_2$ 靠 $q_\varepsilon$ 的近似单位性质（不能用控制收敛）。本篇已在 §11.4 点明。

---

系列：[[10-球与半空间上的Green函数]] · [[09-Green函数的性质]] · [[12-习题课]]
