# 09 Green 函数的性质

> 系列：[[08-Green恒等式与Green表示公式]] · [[10-球与半空间上的Green函数]] · [[11-Poisson积分公式与Dirichlet问题的可解性]]

**定义（Green 函数）**　设 $\Omega\subseteq\mathbb R^n$ 为有界区域，$n\ge2$，边界为 $C^1$。任取 $x\in\Omega$，对 $y\in\overline\Omega\setminus\{x\}$，定义

$$
G(x,y)=\Gamma(x-y)-\gamma(x,y),
$$

其中 $\gamma(x,y)$ 为

$$
\begin{cases}
\Delta_y\gamma(x,y)=0, & y\in\Omega,\\[4pt]
\gamma(x,y)=\Gamma(x-y), & y\in\partial\Omega
\end{cases}
$$

的解，$\Gamma$ 为基本解（见 [[07-基本解与广义函数]]、[[08-Green恒等式与Green表示公式]]）。$G$ 称为区域 $\Omega$ 的 **Green 函数**。

本篇把 $G$ 的四条基本性质逐个证出来。它们既是后面「用 $G$ 解 Dirichlet 问题」的技术基础，也解释了 $G$ 的物理含义。

---

## 9.1 零边值

**性质 1**　对于 $\forall x\in\Omega$、$y\in\partial\Omega$，

$$
G(x,y)=0 .
$$

**证明**　由定义 $G(x,y)=\Gamma(x-y)-\gamma(x,y)$，而在 $\partial\Omega$ 上 $\gamma(x,y)=\Gamma(x-y)$，故 $G=0$。$\square$

---

## 9.2 正性与增长速度估计

**性质 2**　对 $\forall x,y\in\Omega$、$x\ne y$，均有

$$
0<G(x,y)<\begin{cases}
\Gamma(x-y)+\dfrac{1}{2\pi}\log\operatorname{diam}(\Omega), & n=2,\\[8pt]
\Gamma(x-y), & n\ge3,
\end{cases}
$$

且当 $y\to x$ 时，

$$
G(x,y)=\Gamma(x-y)+O(1).
$$

**证明**　分三步。

**第一步：$\gamma$ 有界，故 $G$ 保留了 $\Gamma$ 的奇性。** 由 $\Delta_y\gamma(x,y)=0$，对 $\gamma(\,\cdot\,)$ 在 $\Omega$ 上用极值原理（$y\in\partial\Omega$ 上 $\gamma=\Gamma(x-\cdot)$）：

$$
\min_{y\in\partial\Omega}\Gamma(x-y)\le\gamma(x,y)\le\max_{y\in\partial\Omega}\Gamma(x-y).
$$

也就是说 $\gamma(x,y)$ 在 $y\to x$ 时**有界**，因而

$$
\lim_{y\to x}\frac{G(x,y)}{\Gamma(x-y)}=\lim_{y\to x}\left(1-\frac{\gamma(x,y)}{\Gamma(x-y)}\right)=1 .
$$

**第二步：严格正性（下界）。** 对 $\forall x\in\Omega$，取 $0<r\ll1$ 使 $B_r(x)\Subset\Omega$。则在 $\partial B_r(x)$ 上

$$
\Gamma(x-y)=\Gamma(r)=\begin{cases}
-\dfrac{1}{2\pi}\log r, & n=2,\\[8pt]
\dfrac{1}{(n-2)|S^{n-1}|r^{n-2}}, & n\ge3,
\end{cases}
$$

均有 $\Gamma(r)\gg\max_{y\in\partial\Omega}|\Gamma(x-y)|$（让 $r\to0$ 时左端发散）。于是

$$
G(x,y)\Big|_{\partial B_r(x)}=\Gamma(r)-\gamma(x,y)\ge\Gamma(r)-\max_{y\in\partial\Omega}\Gamma(x-y)>0,
$$

而 $G(x,y)|_{\partial\Omega}=0$。在 $\Omega\setminus\overline{B_r(x)}$ 上 $G(x,\cdot)$ 调和（$\Delta_yG=\Delta_y\Gamma-\Delta_y\gamma=0-0=0$，因为 $y\ne x$）。由极值原理，$G(x,\cdot)$ 的最小值取在边界上，即 $\ge0$；又因为 $G>0$ 于 $\partial B_r(x)$，故 $G$ 不可能恒为常数 $0$，由极值原理的「内部取极值则为常数」可知**内部处处严格为正**：

$$
G(x,y)>0,\qquad\forall y\in\Omega,\ y\ne x .
$$

**第三步：上界。**

- **$n\ge3$ 时**：由第一步，$\gamma(x,y)\ge\min_{y\in\partial\Omega}\Gamma(x-y)>0$（因为 $\Gamma>0$ 处处成立），再用极值原理得 $\gamma>0$ 于 $\Omega$。因而

$$
G(x,y)=\Gamma(x-y)-\gamma(x,y)<\Gamma(x-y).
$$

- **$n=2$ 时**：$\Gamma(x-y)=-\frac{1}{2\pi}\log|x-y|$，由第一步

$$
\gamma(x,y)\ge\min_{y\in\partial\Omega}\Gamma(x-y)=\min_{y\in\partial\Omega}\left(-\frac{1}{2\pi}\log|x-y|\right)
=-\frac{1}{2\pi}\max_{y\in\partial\Omega}\log|x-y|
>-\frac{1}{2\pi}\log\operatorname{diam}(\Omega),
$$

最后一步用 $\max_{y\in\partial\Omega}|x-y|\le\operatorname{diam}(\Omega)$。因而

$$
G(x,y)<\Gamma(x-y)+\frac{1}{2\pi}\log\operatorname{diam}(\Omega).
$$

$\square$

> **补充（第一步与第三步为什么用同一个不等式）**　讲稿把「$\gamma$ 被 $\Gamma$ 在边界上的最大/最小值夹住」用了两次：一次用来断言 $\gamma$ 在 $y\to x$ 处有界（于是在 $x$ 附近 $G\approx\Gamma$），一次用来给出 $G$ 的上界。这个「夹住」就是极值原理的直接推论，而 $\Gamma(x-\cdot)$ 在 $\partial\Omega$ 上的最大最小值都是**显式可算**的——这正是能用它给出定量上下界的原因。

> **补充（第二步中「严格正」要用极值原理的哪一半）**　只由 $\min_{\overline\Omega}G=\min_{\partial\Omega}G$ 只能得到 $G\ge0$（因为 $G|_{\partial\Omega}=0$）。要得到**严格**正，必须再用极值原理的另一半：**非常数的调和函数不能在内部取到最小值**。这里 $G$ 在 $\partial B_r(x)$ 上严格为正，故 $G$ 不可能恒等于 $0$，于是内部不可能取到 $0$。这一点讲稿未点明。

---

## 9.3 对称性

**性质 3**　任取 $x,y\in\Omega$、$x\ne y$，则

$$
G(x,y)=G(y,x).
$$

**证明**　取 $\varepsilon>0$ 足够小，使两闭球 $\overline{B_\varepsilon(x)}$、$\overline{B_\varepsilon(y)}$ 互不相交且均包含于 $\Omega$。在 $\Omega\setminus\bigl(B_\varepsilon(x)\cup B_\varepsilon(y)\bigr)$ 上对 $G(x,z)$ 与 $G(y,z)$ 用 Green 第二公式：

$$
0=\int_{\Omega\setminus(B_\varepsilon(x)\cup B_\varepsilon(y))}\bigl[G(x,z)\Delta_zG(y,z)-G(y,z)\Delta_zG(x,z)\bigr]dz
=\int_{\partial\Omega}\left[\cdots\right]d\sigma_z
-\int_{\partial B_\varepsilon(x)}\left[\cdots\right]d\sigma_z
-\int_{\partial B_\varepsilon(y)}\left[\cdots\right]d\sigma_z
=:I+II+III .
$$

> **观察**　$G(x,z)$、$G(y,z)$ 两个函数均在 $\Omega\setminus\bigl(B_\varepsilon(x)\cup B_\varepsilon(y)\bigr)$ 上调和（因为 $z$ 既不等于 $x$ 也不等于 $y$）。此外，$G(x,z)$ 在 $B_\varepsilon(y)$ 内调和，$G(y,z)$ 在 $B_\varepsilon(x)$ 内调和。

**I（外边界）= 0。** 由

$$
G(x,z)\big|_{z\in\partial\Omega}\equiv G(y,z)\big|_{z\in\partial\Omega}\equiv0,
$$

得 $I=0$。

**II 与 III（内边界）。** 在 $\partial B_\varepsilon(x)$ 上（$\nu_z$ 取小球的单位外法向量，指向球外），

$$
G(x,z)=\Gamma(z-x)-\gamma(x,z)=\Gamma(\varepsilon)-\gamma(x,z),
\frac{\partial G}{\partial\nu_z}(x,z)=\Gamma'(\varepsilon)-\frac{\partial\gamma}{\partial\nu_z}(x,z)=-\frac{1}{|\partial B_\varepsilon|}-\frac{\partial\gamma}{\partial\nu_z}(x,z).
$$

于是

$$
II=-\int_{\partial B_\varepsilon(x)}\bigl(\Gamma(\varepsilon)-\gamma(x,z)\bigr)\frac{\partial G}{\partial\nu_z}(y,z)\,d\sigma_z
+\int_{\partial B_\varepsilon(x)}G(y,z)\left(-\frac{1}{|\partial B_\varepsilon|}-\frac{\partial\gamma}{\partial\nu_z}(x,z)\right)d\sigma_z .
$$

其中三项分别处理：

- $\Gamma(\varepsilon)\displaystyle\int_{\partial B_\varepsilon(x)}\frac{\partial G}{\partial\nu_z}(y,z)\,d\sigma_z=0$，因为 $G(y,\cdot)$ 在 $B_\varepsilon(x)$ 内调和，故其法向导数在球面上的积分为 $\displaystyle\int_{B_\varepsilon(x)}\Delta_zG(y,z)dz=0$；
- 含 $\gamma$ 的左右两项合起来是

$$
\int_{\partial B_\varepsilon(x)}\left[\gamma(x,z)\frac{\partial G(y,z)}{\partial\nu_z}-G(y,z)\frac{\partial\gamma(x,z)}{\partial\nu_z}\right]d\sigma_z
=\int_{B_\varepsilon(x)}\bigl[\gamma\,\Delta_zG(y,\cdot)-G(y,\cdot)\,\Delta_z\gamma\bigr]dz=0
$$

> （Green 第二公式 + 两者在 $B_\varepsilon(x)$ 内调和）；

- 剩下的就是

$$
II=-\frac{1}{|\partial B_\varepsilon|}\int_{\partial B_\varepsilon(x)}G(y,z)\,d\sigma_z .
$$

（这正是「小球面平均」，取 $\varepsilon\to0$ 时会给出 $-G(y,x)$。）类似地，

$$
III=\frac{1}{|\partial B_\varepsilon|}\int_{\partial B_\varepsilon(y)}G(x,z)\,d\sigma_z .
$$

**取极限。** 令 $\varepsilon\to0$，由 $G$ 的连续性（奇性为 $\Gamma\sim O(|z-x|^{2-n})$，其在球面上的平均收敛到函数值）

$$
0=I+II+III\ \longrightarrow\ 0-G(y,x)+G(x,y),
$$

故 $G(x,y)=G(y,x)$。$\square$

> **要点**　对称性证明的机制是「**用 Green 第二公式把一个反对称的量夹在两个小球的积分之间**」。注意内边界项 $II$、$III$ 的符号相反（因为两球的外法向量方向相反），最终给出 $G(x,y)-G(y,x)=0$。物理上这就是**互易性**（reciprocity）：在 $x$ 放电荷、在 $y$ 测电势，与反过来做，结果相同。

---

## 9.4 边界通量

**性质 4**　对 $\forall x\in\Omega$，均成立

$$
\int_{\partial\Omega}\frac{\partial G(x,y)}{\partial\nu_y}\,d\sigma_y=-1 .
$$

**证明**　取 $\varepsilon>0$，$\text{s.t. } B_\varepsilon(x)\Subset\Omega$。在 $\Omega\setminus B_\varepsilon(x)$ 上 $G(x,\cdot)$ 调和，由散度定理

$$
0=\int_{\Omega\setminus B_\varepsilon(x)}\Delta_yG(x,y)\,dy
=\int_{\partial\Omega}\frac{\partial G}{\partial\nu_y}(x,y)\,d\sigma_y
-\int_{\partial B_\varepsilon(x)}\frac{\partial G}{\partial\nu_y}(x,y)\,d\sigma_y .
$$

> 第二项前的负号来自：作为 $\Omega\setminus B_\varepsilon(x)$ 的外法向量，在 $\partial B_\varepsilon(x)$ 上应取小球的 $- \nu_B$；此处已把 $\nu_y$ 统一记为小球的单位外法向量 $\nu_B=\dfrac{y-x}{|y-x|}$，故符号写进减号里。

而 $y\in\partial B_\varepsilon(x)$ 时，

$$
\frac{\partial G}{\partial\nu_y}(x,y)=\frac{\partial\Gamma(x-y)}{\partial\nu_y}-\frac{\partial\gamma(x,y)}{\partial\nu_y}
=-\frac{1}{|\partial B_\varepsilon(x)|}-\frac{\partial\gamma(x,y)}{\partial\nu_y},
$$

于是

$$
-\int_{\partial B_\varepsilon(x)}\frac{\partial G}{\partial\nu_y}(x,y)\,d\sigma_y
=1+\int_{\partial B_\varepsilon(x)}\frac{\partial\gamma}{\partial\nu_y}\,d\sigma_y
=1+\int_{B_\varepsilon(x)}\Delta_y\gamma(x,y)\,dy=1 .
$$

代回得

$$
0=\int_{\partial\Omega}\frac{\partial G(x,y)}{\partial\nu_y}\,d\sigma_y-1
\quad\Longrightarrow\quad
\int_{\partial\Omega}\frac{\partial G(x,y)}{\partial\nu_y}\,d\sigma_y=-1 .
$$

$\square$

> **要点**　性质 4 说明 $-\dfrac{\partial G(x,y)}{\partial\nu_y}$ 在 $\partial\Omega$ 上是一个**归一化的核**（总积分为 $1$）。这一点在 [[11-Poisson积分公式与Dirichlet问题的可解性]] 里至关重要：验证 $u(x)=-\displaystyle\int_{\partial\Omega}g(y)\dfrac{\partial G}{\partial\nu_y}d\sigma_y$ 在边界上取到 $g$ 时，正是靠「$\int(-\partial G/\partial\nu_y)=1$」把 $g(x_0)$ 写成积分，从而把 $u(x)-g(x_0)$ 化成 $g(y)-g(x_0)$ 的积分。

---

## 9.5 静电学模型

设 $\Omega\subseteq\mathbb R^3$ 为一有界区域，$\Omega$ 中有电荷密度函数 $\rho(x)\ge0$，则由 Maxwell 方程组

$$
\begin{cases}
\operatorname{div}E=\dfrac{\rho}{\varepsilon_0}, & \varepsilon_0>0\text{ 为真空介电常数},\\[8pt]
\operatorname{curl}E=0, &
\end{cases}
$$

（当 $\Omega$ 单连通时）存在 $u$ 使得

$$
-\Delta u=\frac{\rho}{\varepsilon_0},\qquad x\in\Omega .
$$

**于是 Green 函数具有静电学意义**：设 $\partial\Omega$ 为接地导体，在 $x\in\Omega$ 处放置电量为 $q$ 的电荷，则电势为

$$
u(y)=\frac{q}{\varepsilon_0}\cdot G(x,y).
$$

> **补充（四条性质的物理读法）**
>
> | 性质 | 数学表述 | 物理含义 |
> |---|---|---|
> | 零边值 | $G(x,y)=0$ 于 $y\in\partial\Omega$ | 边界是接地的理想导体，电势为零 |
> | 正性 | $G>0$ | 正电荷在接地腔内的电势处处为正 |
> | 奇性 | $G=\Gamma(x-y)+O(1)$ | 靠近点电荷时，$G$ 与自由空间基本解同奇性（单位点电荷的场） |
> | 对称性 | $G(x,y)=G(y,x)$ | 互易性：源与观测点互换，结果不变 |
> | 边界通量 | $\int_{\partial\Omega}\partial G/\partial\nu_y=-1$ | 穿过整个导体边界的电通量归一 |

> **要点**　Green 函数的「Green」与 Green 公式同源，但它的真正内容是「**把边界条件吸收进核里**」：$\Gamma(x-y)$ 是自由空间的响应，$\gamma(x,y)$ 是边界（接地导体）产生的**感应场**，$G=\Gamma-\gamma$ 是两者之差，正好让边界上的电势归零。

---

## 订正

1. **性质 2 的正性证明只「用极值原理」，不足以得到严格不等号。** 由 $\min_{\overline\Omega}G=\min_{\partial\Omega}G$ 只能得 $G\ge0$（因为 $G|_{\partial\Omega}=0$）。要得 $G>0$，必须再用极值原理的「非常数调和函数不能在内部取到最小值」。讲稿原文为「对 $G(x,y)$ 在 $\Omega\setminus B_r(x)$ 上使用极值原理知 $G(x,y)>0$」，跳过了这一半。本篇已在**补充**中补出。
2. **性质 2 中 $n=2$ 的上界，讲稿省略了 $\max_{y\in\partial\Omega}|x-y|\le\operatorname{diam}(\Omega)$ 这一步。** 正是靠它才能把 $\min_{y\in\partial\Omega}\bigl(-\frac{1}{2\pi}\log|x-y|\bigr)$ 换成 $-\frac{1}{2\pi}\log\operatorname{diam}(\Omega)$。本篇已补。
3. **性质 3 的化简一步跨得较大。** 讲稿把 $II$ 的三项直接合并为 $II=-\dfrac{1}{|\partial B_\varepsilon|}\displaystyle\int_{\partial B_\varepsilon(x)}G(y,z)d\sigma_z$，未说明前两项为何消失（$\Gamma(\varepsilon)$ 的通量积分为零、含 $\gamma$ 的两项由 Green 第二公式相消）。本篇已逐项补出。
4. **性质 4 中 $\nu_y$ 的取法前后不一致。** 讲稿先说「在 $\partial B_\varepsilon(x)$ 上」再直接把 $\partial G/\partial\nu_y$ 写成 $-1/|\partial B_\varepsilon(x)|-\partial\gamma/\partial\nu_y$，其间隐含地把 $\nu_y$ 取成了**小球的外法向量**，而散度定理给出的是 $- \nu_B$。本篇已在括号中写明约定，避免符号来源不明。
5. **性质 2 的假设中「$y\ne x$」不可省。** 讲稿写成 $\forall x,y\in\Omega$，紧跟着写 $x\ne y$，格式上像是附加条件；实际上 $G$ 在 $y=x$ 处没有定义（$\Gamma$ 发散），$y\ne x$ 是定义域的限制而非额外假设。本篇已调整表述。

---

系列：[[08-Green恒等式与Green表示公式]] · [[10-球与半空间上的Green函数]] · [[11-Poisson积分公式与Dirichlet问题的可解性]]
