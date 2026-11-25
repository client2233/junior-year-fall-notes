# 01 从万有引力到 Poisson 方程

> 系列：[[02-边值问题的开始·三类边界条件与Dirichlet原理]] · [[03-调和函数·平均值性质与极值原理]] · [[07-基本解与广义函数]]

本课程研究的是**有具体物理背景的偏微分方程**。回忆常微分方程：设 $x(t)$ 是随时间变化的未知量，则

$$
\dot x(t)=f\bigl(x(t),t\bigr)
$$

只含一个自变量 $t$，未知函数 $x$ 的导数只有对 $t$ 的。若要描述一块薄板的温度，同一时刻不同位置的温度也会不同，记 $u(x,y,t)$ 为位置 $(x,y)$ 在时刻 $t$ 的温度，则温度随时间与随位置的变化分别由对 $t$ 与对 $x,y$ 的偏导数描述。含未知函数偏导数的方程称为**偏微分方程**，热方程就是一个典型例子。

本篇从万有引力出发，导出本课程的两条主线方程：Poisson 方程与 Laplace 方程。

---

## 1.1 牛顿运动定律与引力场

牛顿第二定律把力定义为动量的瞬时变化率：

$$
F=\frac{d(mv)}{dt}.
$$

记 $a=\dfrac{dv}{dt}$ 为加速度。质量 $m$ 通常恒定，此时 $F=ma$，受力大小恰好对应加速度。

**万有引力**　设在点 $y\in\mathbb R^3$ 处有质量 $\rho(y)$，则在点 $x\ne y$ 处、质量为 $m(x)$ 的质点受到的引力为

$$
F_y(x)=-m(x)\,\frac{G\rho(y)}{|x-y|^3}\,(x-y).\tag{1.1}
$$

这个公式里有四点值得拆开看：

1. 除去试验质量 $m(x)$，可考虑

$$
g_y(x)=-\frac{G\rho(y)}{|x-y|^3}\,(x-y).
$$

这是点 $y$ 处的质量 $\rho(y)$ 在 $\mathbb R^3\setminus\{y\}$ 上产生的向量值函数，称为**引力加速度场**。有了它，在不同点、不同质量对应的引力相乘即可得出。

2. 若假设质量集中在点 $y$、在 $y$ 之外没有质量，则可以把 $\rho(y)$ 看作空间中的函数，称为**质量密度函数**。

3. $\dfrac{x-y}{|x-y|}$ 是从 $y$ 指向 $x$ 的单位方向向量，故 (1.1) 中的负号表示引力**指向** $y$。

4. $G>0$ 是源于正比例关系的常数，称为**万有引力常数**。

### 问题 1：为什么是平方反比

固定彼此质量，为什么

$$
|F_y(x)|\propto\frac{1}{|x-y|^2}\;?
$$

考虑以下两个假设：

- **假设 1（球对称）**：质点对周围空间的作用是球对称的。把 $\rho(y)$ 看作质点，直观上没道理不同方向产生的力不一样。
- **假设 2（通量均匀）**：质点对周围空间所产生的作用是均匀的，即穿过不同半径球面的引力场通量相同（两球面之间没有其他质量源）。

以 $y$ 为球心，比较半径 $r$ 与 $R$ 的球面。三维空间中半径 $R$ 的球面面积为 $4\pi R^2$，通量相同而面积按 $R^2$ 增长，故

$$
|g_y(x)|\propto\frac{1}{4\pi R^2},\qquad |x-y|=R.
$$

对固定的试验质量，$|F_y(x)|$ 也具有相同的平方反比规律。

类比地，在相应的对称性与通量假设下：二维空间中引力大小正比于 $1/(2\pi R)$，四维空间中正比于 $1/R^3$。

> **要点**　这里的「次数」$n-1$ 恰好是 $\mathbb R^n$ 中球面面积的增长阶。它正是 Laplace 方程基本解 $\Gamma(x)\propto|x-y|^{2-n}$ 中出现的次数——两组指数互补，后面 [[07-基本解与广义函数]] 会看到这不是巧合。

### 问题 2：质量体的引力场

实际计算引力场时需要考虑一个质量体 $V$。设 $\rho$ 为质量密度，并在 $V$ 外延拓为零：

$$
\rho(x)=\begin{cases}\text{质量密度},& x\in V,\\[2pt] 0,& x\notin V.\end{cases}
$$

最朴素的想法是把各处产生的引力迭加起来：

$$
g_V(x)=\int_V g_y(x)\,dy=-G\int_{\mathbb R^3}\frac{\rho(y)}{|x-y|^3}\,(x-y)\,dy.\tag{1.2}
$$

其中 $g_y(x)\,dy$ 表示体积元 $dy$ 中的质量 $\rho(y)\,dy$ 产生的引力加速度。这个积分通常很难直接计算。

**关键观察**　当 $x\ne y$ 时，

$$
\nabla_x\frac{1}{|x-y|}=-\frac{x-y}{|x-y|^3},
$$

即

$$
g_y(x)=\nabla_x\!\left(\frac{G\rho(y)}{|x-y|}\right).
$$

于是若记

$$
u(x)=\int_{\mathbb R^3}\frac{G\rho(y)}{|x-y|}\,dy,\qquad g_V(x)=\nabla u(x).\tag{1.3}
$$

这种形式的积分所定义的函数称为**卷积函数**。为使用下文的经典微分与积分公式，暂设密度 $\rho\in C_c^\infty(\mathbb R^3)$。

> **补充（关键观察的验证）**　$|x-y|=\Bigl(\sum_i(x_i-y_i)^2\Bigr)^{1/2}$，故 $\partial_{x_i}|x-y|=\dfrac{x_i-y_i}{|x-y|}$，即 $\nabla_x|x-y|=\dfrac{x-y}{|x-y|}$。再对 $|x-y|^{-1}$ 用链式法则：

$$
\nabla_x\frac{1}{|x-y|}=-\frac{1}{|x-y|^2}\nabla_x|x-y|=-\frac{x-y}{|x-y|^3}.
$$

---

## 1.2 势函数的物理学含义

先看点质量。记

$$
u_y(z)=\frac{G\rho(y)}{|z-y|}.
$$

取从 $x$ 经过 $p$ 到 $y$ 的曲线 $\gamma$，使

$$
\gamma(0)=x,\qquad \gamma(1)=p,\qquad \gamma(2)=y,
$$

只考虑避开源点 $y$ 的曲线段 $\gamma|_{[0,1]}$。由于 $g_y=\nabla u_y$，引力对单位质量所做的功为

$$
\int_0^1\bigl\langle g_y(\gamma(t)),\gamma'(t)\bigr\rangle\,dt
=\int_0^1\frac{d}{dt}u_y(\gamma(t))\,dt
=u_y(p)-u_y(x)
=\frac{G\rho(y)}{|p-y|}-\frac{G\rho(y)}{|x-y|}.
$$

因此 $-u_y$ 是点 $y$ 处质量对单位质量的引力势能；对一般质量体，单位质量的引力势能为 $-u$。

类比地，若把相应维数中径向引力场的比例常数仍记为 $G$，则二维空间中可写为

$$
g_V(x)=\nabla_x\left(-G\int_{\mathbb R^2}\rho(y)\log|x-y|\,dy\right),
$$

四维空间中可写为

$$
g_V(x)=\nabla_x\left(\frac{G}{2}\int_{\mathbb R^4}\frac{\rho(y)}{|x-y|^2}\,dy\right).
$$

更简单地记 $g=\nabla u$。这个看法并没有使直接计算变得容易，但可以迂回讨论。

> **要点**　注意二维与四维的符号差异：$n\ge3$ 时基本解是 $+\dfrac{1}{(n-2)|S^{n-1}|\,|x-y|^{n-2}}$，而 $n=2$ 时是 $-\dfrac{1}{2\pi}\log|x-y|$。二维之所以多一个负号，是因为 $\log r$ 与 $r^{2-n}$ 的单调性相反。

---

## 1.3 引力场的通量与 Poisson 方程

**定理 1.1**　沿用前面的假设，设 $\rho\in C_c^\infty(\mathbb R^3)$ 且 $\rho\ge0$，令

$$
u(x):=\int_{\mathbb R^3}\frac{G\rho(y)}{|x-y|}\,dy.
$$

则

$$
-\Delta u=4\pi G\rho.\tag{1.4}
$$

**证明**　我们证明：对任意具有 $C^1$ 边界的有界区域 $\Omega\subset\mathbb R^3$，均有

$$
\int_\Omega\bigl(\Delta u+4\pi G\rho\bigr)\,dx=0.
$$

由区域 $\Omega$ 的任意性与被积函数的连续性即得结论。

**第一步：把体积分化成通量。** 记 $\nu(x)$ 为单位外法向量，$g=\nabla u$。一方面，由 Gauss–Green 公式（见 1.3.1），

$$
\int_{\partial\Omega}\langle g,\nu\rangle\,d\sigma_x=\int_\Omega\operatorname{div}g\,dx=\int_\Omega\Delta u\,dx.\tag{1.5}
$$

**第二步：换一种算法算同一个通量。** 利用 (1.3) 并交换积分次序，

$$
\int_{\partial\Omega}\langle\nabla u,\nu\rangle\,d\sigma_x
=\int_{\partial\Omega}\left\langle\int_{\mathbb R^3}G\rho(y)\nabla_x\frac{1}{|x-y|}\,dy,\;\nu(x)\right\rangle d\sigma_x
=\int_{\partial\Omega}\!\int_{\mathbb R^3}G\rho(y)\left\langle\nabla_x\frac{1}{|x-y|},\nu(x)\right\rangle dy\,d\sigma_x.
$$

把对 $y$ 的积分按 $y\in\Omega$ 与 $y\in\mathbb R^3\setminus\Omega$ 拆成两段（$\partial\Omega$ 的三维体积为零，不影响对 $y$ 的积分）：

$$
=I+II.
$$

**第三步：算 $I$。** 固定 $y\in\Omega$，取 $\varepsilon>0$ 使 $B_\varepsilon(y)\subset\Omega$。记 $\nu_B$ 为小球 $\partial B_\varepsilon(y)$ 的单位外法向量。在挖孔域 $\Omega\setminus B_\varepsilon(y)$ 上用 Gauss–Green 公式，其边界由 $\partial\Omega$ 与 $\partial B_\varepsilon(y)$ 组成，而在后一片上区域的外法向量恰为 $-\nu_B$，故

$$
\int_{\partial\Omega}\left\langle\nabla_x\frac{1}{|x-y|},\nu(x)\right\rangle d\sigma_x
=\int_{\Omega\setminus B_\varepsilon(y)}\operatorname{div}_x\nabla_x\frac{1}{|x-y|}\,dx
+\int_{\partial B_\varepsilon(y)}\left\langle\nabla_x\frac{1}{|x-y|},\nu_B\right\rangle d\sigma_x.
$$

当 $x\ne y$ 时 $\Delta_x\dfrac{1}{|x-y|}=0$，故体积分为零。又因在 $\partial B_\varepsilon(y)$ 上有 $\nu_B=\dfrac{x-y}{|x-y|}$，于是被积函数是常数：

$$
\left\langle-\frac{x-y}{|x-y|^3},\;\nu_B\right\rangle
=-\frac{|x-y|^2}{|x-y|^4}\Bigg|_{|x-y|=\varepsilon}=-\frac{1}{\varepsilon^2},
$$

从而

$$
\int_{\partial B_\varepsilon(y)}\left\langle\nabla_x\frac{1}{|x-y|},\nu_B\right\rangle d\sigma_x=-\frac{1}{\varepsilon^2}\cdot4\pi\varepsilon^2=-4\pi.
$$

代回并交换积分次序得

$$
I=\int_\Omega G\rho(y)\cdot(-4\pi)\,dy=-4\pi G\int_\Omega\rho(y)\,dy.
$$

**第四步：算 $II$。** 此时 $y\notin\Omega$，核函数 $\dfrac{1}{|x-y|}$ 作为 $x$ 的函数在 $\overline\Omega$ 上光滑且**调和**，故

$$
\int_{\partial\Omega}\left\langle\nabla_x\frac{1}{|x-y|},\nu(x)\right\rangle d\sigma_x=\int_\Omega\Delta_x\frac{1}{|x-y|}\,dx=0,
$$

于是 $II=0$。

**第五步：合并。** 由 (1.5) 与上面两段，

$$
\int_\Omega\Delta u\,dx=-4\pi G\int_\Omega\rho(y)\,dy,
$$

即

$$
\int_\Omega\bigl(\Delta u+4\pi G\rho\bigr)\,dx=0.
$$

由 $\Omega$ 的任意性，被积函数恒为零，即 $-\Delta u=4\pi G\rho$。$\square$

**Poisson 方程与 Laplace 方程**　(1.4) 称为 **Poisson 方程**。若在所考察区域内 $\rho\equiv0$，则得到

$$
\Delta u=0,
$$

称为 **Laplace 方程**，也称**调和方程**。

> **要点**　证明的实质是「把 $\Delta$ 从待求的 $u$ 挪到显式已知的核 $1/|x-y|$ 上」。核在 $x=y$ 处有奇性，所以必须挖孔、单独处理小球面上的通量——这个 $-4\pi$ 正是 Poisson 方程里 $4\pi$ 的来源。凡是要证明「某个卷积核确实是基本解」，用的都是同一套挖孔手法。

### 1.3.1 Gauss–Green 公式

设 $\Omega\subset\mathbb R^3$ 为具有 $C^1$ 边界的有界区域，$F=(F_1,F_2,F_3)\in C^1(\overline\Omega;\mathbb R^3)$ 为向量场，$\nu$ 为单位外法向量。则

$$
\int_\Omega\operatorname{div}F\,dx=\int_{\partial\Omega}\langle F,\nu\rangle\,d\sigma_x,\qquad
\operatorname{div}F=\sum_{i=1}^{3}\frac{\partial F_i}{\partial x_i}.\tag{1.7}
$$

**证明概要**　对区域 $\Omega$ 作分割：

1. 先在长方体 $I=[a,b]\times[c,d]\times[e,f]$ 上证明；
2. 再在与边界相交的块上证明；
3. 最后把每个块上的结果相加。

> **补充（长方体情形的完整证明）**　只需对单个分量 $\partial_iF_i$ 证明。以 $\partial_3F_3$ 为例，记 $\Omega=I$，

$$
\int_I\frac{\partial F_3}{\partial x_3}\,dx
=\int_{[a,b]\times[c,d]}\left(\int_e^f\frac{\partial F_3}{\partial x_3}\,dx_3\right)dx_1dx_2
=\int_{[a,b]\times[c,d]}\bigl(F_3(\cdot,\cdot,f)-F_3(\cdot,\cdot,e)\bigr)\,dx_1dx_2.
$$

在顶面 $x_3=f$ 上外法向量是 $\nu=(0,0,1)$，$\langle F,\nu\rangle=F_3$；在底面 $x_3=e$ 上是 $\nu=(0,0,-1)$，$\langle F,\nu\rangle=-F_3$；四个侧面 $\nu_3=0$，无贡献。于是上式恰为

$$
\int_{\partial I}\langle F,\nu\rangle\,d\sigma_x=\int_{\partial I}F_3\nu_3\,d\sigma_x .
$$

对 $i=1,2$ 同理，相加即得长方体上的 (1.7)。一般 $C^1$ 区域用单位分解把被积函数局部化，每片上的坐标图把区域拉直为长方体，再用上面结论并求和——这正是上面三步的作用。

---

## 1.4 场的数学模型总结

万有引力的加速度场是

$$
g=\nabla u,\qquad u(x)=\int_{\mathbb R^3}\frac{G\rho(y)}{|x-y|}\,dy,
$$

即卷积形式的表达。相应的**全空间模型**为

$$
\begin{cases}
-\Delta u=4\pi G\rho, & \text{在 }\mathbb R^3\text{ 中},\\[4pt]
\lim\limits_{|x|\to\infty}u(x)=0, & \text{(无穷远条件)}.
\end{cases}\tag{1.8}
$$

**无穷远处极限的验证**　设质量体 $V\subset B_R(0)$，且 $\rho$ 在 $V$ 外为零。当 $|x|>R$ 时，由 $|x-y|\ge|x|-|y|\ge|x|-R$ 得

$$
0\le u(x)=\int_{\mathbb R^3}\frac{G\rho(y)}{|x-y|}\,dy=\int_V\frac{G\rho(y)}{|x-y|}\,dy
\le\frac{1}{|x|-R}\int_VG\rho(y)\,dy\longrightarrow0\qquad(|x|\to\infty).
$$

**类似的例子：固定电荷的静电场**　单连通区域中固定的电荷密度函数 $\rho$ 定义的静电场 $E$。取静电比例常数为 $1$，则满足

$$
\operatorname{curl}E=0,\qquad \operatorname{div}E=\rho.
$$

旋度为 $0$ 意味着 $E$ 是梯度场（可通过路径积分定义出原函数），从而

$$
E=-\nabla u,\qquad -\Delta u=\rho.
$$

> **要点**　同一套数学结构（$-\Delta u=$ 源）在不同的物理外衣下反复出现。注意符号约定：引力场 $g=+\nabla u$，而静电场 $E=-\nabla u$，这是后面 [[02-边值问题的开始·三类边界条件与Dirichlet原理]] 中讨论唯一性与边值条件时的固定约定。

---

## 订正

本篇（对应讲稿第一章 1.1–1.4）经逐式复核，未发现排印错误。两点说明：

1. §1.3.1 的 Gauss–Green 公式，讲稿明确只保留了证明步骤（「边界块的处理及求和过程尚未展开」），本篇已在**补充**中补全长方体情形并说明一般区域的处理路径。
2. (1.2) 中 $g_y(x)$ 的 $\rho(y)$ 与质量体情形下的「密度函数」同符号，是讲稿的有意记号复用（讲稿 1.1 第 2 条已说明「可以把 $\rho(y)$ 看作空间中的函数」），初学者易误读为点质量密度，此处不作订正。

---

系列：[[02-边值问题的开始·三类边界条件与Dirichlet原理]] · [[03-调和函数·平均值性质与极值原理]] · [[07-基本解与广义函数]]
