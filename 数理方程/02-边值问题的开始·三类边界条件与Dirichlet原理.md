# 02 边值问题的开始·三类边界条件与 Dirichlet 原理

> 系列：[[01-从万有引力到Poisson方程]] · [[03-调和函数·平均值性质与极值原理]] · [[04-调和函数·Harnack不等式与内部梯度估计]]

上一章从万有引力导出了 Poisson 方程与 Laplace 方程。本篇先看另外两个产生这类方程的例子，再讨论一个核心问题：**仅仅给出区域内的方程，能否确定所求的场或势函数？** 结论是不能——必须补上边界信息，由此引出三类边值问题。最后给出「解方程」的第一个非构造性工具：变分原理。

---

## 2.1 另外两个场的例子

### 2.1.1 固定电荷的静电场

设 $\Omega\subseteq\mathbb R^3$ 为单连通区域，区域内有固定的电荷密度函数 $\rho(x_1,x_2,x_3)$。取静电比例常数为 $1$，则静电场 $E$ 满足

$$
\operatorname{div}E=\rho,\qquad \operatorname{curl}E=0.
$$

由 $\Omega$ 的单连通性与 $\operatorname{curl}E=0$，可以引入**静电势** $u$，使

$$
E=-\nabla u.
$$

代入散度方程，得到

$$
-\Delta u=\rho.\tag{2.1}
$$

这里 $u$ 表示静电势，采用 $E=-\nabla u$ 的约定；而上一章的正牛顿势满足 $g=\nabla u$。两者相差一个整体符号。

> **要点**　单连通性不可省。在有多连通区域上，无旋场只在「每个闭合回路的环量为零」时才是有势场；一般情形下 $E$ 只能写成 $-\nabla u$ 加一个无散无旋的调和场。这也是后面 [[12-习题课]] 第 1 题要用可缩性（比单连通更强）去构造势函数的原因。

### 2.1.2 不可压缩流体的速度场

设 $\Omega\subseteq\mathbb R^3$ 为单连通区域，流体的速度场为 $v(x_1,x_2,x_3)$。若速度场可以写成

$$
v=-\nabla\varphi,
$$

且流体不可压缩，即

$$
\operatorname{div}v=0,
$$

则

$$
\Delta\varphi=0.
$$

因此，在这一**势流**假设下，速度势 $\varphi$ 满足 Laplace 方程。

> **要点**　三个来源给出同一个方程骨架：引力与电场给 Poisson 方程（右端是「源」），不可压缩势流给 Laplace 方程（区域内无源）。所以整门课的重心是：在一个有界区域内解 $-\Delta u=f$，并在边界上补足条件。

---

## 2.2 从万有引力看解的唯一性

以万有引力模型为例。即使规定了 $\rho$，**势函数对应的方程也不一定有唯一解**。设 $\Omega\subset\mathbb R^3$ 且 $\Omega\ne\mathbb R^3$，只给定 $\Omega$ 内的 $\rho$，则

$$
-\Delta u=4\pi G\rho\qquad\text{在 }\Omega\text{ 中}
$$

的解不唯一：若 $u_1$ 与 $u_2$ 均满足这个方程，由求导的线性性，

$$
-\Delta(u_1-u_2)=0,
$$

即 $u_1$ 与 $u_2$ 相差一个**调和函数**。

**例**　$x_1^2-x_2^2$ 就是一个调和函数：

$$
\Delta(x_1^2-x_2^2)=\frac{\partial^2}{\partial x_1^2}(x_1^2)-\frac{\partial^2}{\partial x_2^2}(x_2^2)=2-2=0.
$$

**思考**　如何确定唯一的解？

从物理角度看，需要给出边界上的信息。下面设 $\Omega$ 具有光滑边界，$\nu$ 为单位外法向量。

> **要点**　「相差一个调和函数」这件事把两个问题合成了一个：**Poisson 方程的唯一性完全归结为调和函数的唯一性**。所以第 3–6 篇花很大力气讨论调和函数（平均值性质、极值原理、Harnack 不等式、Liouville 定理），正是为了后面能一言断定「相差调和函数」的那部分一定是零。

---

## 2.3 Dirichlet、Neumann 与 Robin 边值问题

### 2.3.1 Dirichlet 边界条件

固定边界上的势函数**值**，考虑

$$
\begin{cases}
-\Delta u=4\pi G\rho, & \text{在 }\Omega\text{ 中},\\[4pt]
u=\Phi, & \text{在 }\partial\Omega\text{ 上}.
\end{cases}\tag{2.2}
$$

这称为 **Dirichlet 边值问题**，也称**第一边值问题**。

### 2.3.2 Neumann 边界条件

固定边界上的向外**通量**，考虑

$$
\begin{cases}
-\Delta u=4\pi G\rho, & \text{在 }\Omega\text{ 中},\\[4pt]
\dfrac{\partial u}{\partial\nu}=g, & \text{在 }\partial\Omega\text{ 上}.
\end{cases}\tag{2.3}
$$

这称为 **Neumann 边值问题**，也称**第二边值问题**。这里 $g$ 表示给定的标量边界数据，且

$$
\frac{\partial u}{\partial\nu}=\langle\nabla u,\nu\rangle.
$$

> **注**　Neumann 条件仍不能确定势函数的加法常数：若 $u$ 是解，则 $u+C$ 满足同样的方程和边界条件。

### 2.3.3 Robin 边界条件

还有第三种边值问题——**Robin 边值问题**：

$$
\begin{cases}
-\Delta u=4\pi G\rho, & \text{在 }\Omega\text{ 中},\\[4pt]
\dfrac{\partial u}{\partial\nu}+cu=g, & \text{在 }\partial\Omega\text{ 上}.
\end{cases}\tag{2.4}
$$

这里取 $c>0$。**若不限制 $c$ 的符号，唯一性可能不成立。** 例如在单位球 $B_1\subset\mathbb R^n$ 上，

$$
\begin{cases}
-\Delta u=0, & \text{在 }B_1\text{ 中},\\[4pt]
\dfrac{\partial u}{\partial\nu}-u=0, & \text{在 }\partial B_1\text{ 上}
\end{cases}
$$

有 $0,x_1,x_2,\dots,x_n$ 这些不同的解。

> **补充（反例的验证）**　取 $u(x)=x_i$。首先 $\Delta x_i=0$ 于 $B_1$ 中。其次在 $\partial B_1$ 上 $\nu(x)=x$（单位外法向量），故

$$
\frac{\partial u}{\partial\nu}=\langle\nabla x_i,\nu\rangle=\langle e_i,x\rangle=x_i=u,
$$

即 $\dfrac{\partial u}{\partial\nu}-u=0$。所以 $c=-1$ 时齐次 Robin 问题有无穷多解，正是「不限制 $c$ 的符号」的后果。第三边值问题源于热方程（边界与外界的热交换）。

### 2.3.4 三类边值问题对照

| 问题 | 区域内方程 | 边界条件 | 唯一性 |
|---|---|---|---|
| Dirichlet（第一边值） | $-\Delta u=4\pi G\rho$ | $u=\Phi$ 于 $\partial\Omega$ | 若存在则唯一（[[03-调和函数·平均值性质与极值原理]]） |
| Neumann（第二边值） | 同上 | $\dfrac{\partial u}{\partial\nu}=g$ 于 $\partial\Omega$ | 若存在则至多相差一个常数（[[05-下调和函数·弱极值原理与Hopf边界点引理]]） |
| Robin（第三边值） | 同上 | $\dfrac{\partial u}{\partial\nu}+cu=g$，$c>0$ | 若存在则唯一；$c$ 变号时不成立 |

---

## 2.4 变分法与 Dirichlet 原理

设 $\Omega\subset\mathbb R^n$ 为具有光滑边界的有界区域，考虑函数的集合

$$
U:=\bigl\{u\in C^2(\Omega)\cap C^1(\overline\Omega)\ :\ u|_{\partial\Omega}=g\bigr\},
$$

以及它上面的**能量泛函**（也称 **Dirichlet 能量**）

$$
E(u):=\frac12\int_\Omega|\nabla u|^2\,dx.
$$

$E(u)$ 的物理含义：在引力场模型中 $E(u)=\dfrac12\displaystyle\int_\Omega|g|^2dx_1dx_2dx_3$ 是场能；若 $\Omega\subset\mathbb R^2$ 是一张鼓面、$u(x_1,x_2)$ 是高度函数，则 $E(u)$ 是弹性势能。下面仅讨论

$$
\begin{cases}
-\Delta u=0, & \text{在 }\Omega\text{ 中},\\[4pt]
u=g, & \text{在 }\partial\Omega\text{ 上}
\end{cases}\qquad(*)
$$

在 $U$ 中的解。

**定理 2.1（Dirichlet 原理）**　$u\in U$ 为 $E$ 在 $U$ 中的极小值点，当且仅当

$$
\begin{cases}
-\Delta u=0, & \text{在 }\Omega\text{ 中},\\[4pt]
u=g, & \text{在 }\partial\Omega\text{ 上}.
\end{cases}
$$

**证明**　记

$$
U_0:=\bigl\{w\in C^2(\Omega)\cap C^1(\overline\Omega)\ :\ w|_{\partial\Omega}=0\bigr\},
$$

则 $U_0$ 是线性空间，且对任意 $u\in U$ 有 $u+U_0=U$（即 $u-v\in U_0$）。

**（$\Rightarrow$）设 $u\in U$ 为极小值点，证明 $\Delta u=0$。** 任取 $w\in U_0$，作扰动 $u_t:=u+tw$（仍属于 $U$）。计算

$$
\frac{d}{dt}E(u+tw)=\frac{d}{dt}\frac12\int_\Omega|\nabla(u+tw)|^2dx
=\int_\Omega\langle\nabla(u+tw),\nabla w\rangle\,dx.
$$

因为 $u$ 为极小值点，所以

$$
\left.\frac{d}{dt}E(u+tw)\right|_{t=0}=0,
$$

也即

$$
\int_\Omega\langle\nabla u,\nabla w\rangle\,dx=0\qquad\text{对任意 }w\in U_0.\tag{*}
$$

**先取紧支撑于 $\Omega$ 的 $w\in U_0$**，分部积分得

$$
0=\int_\Omega\langle\nabla u,\nabla w\rangle\,dx
=\int_\Omega\left(\operatorname{div}(w\nabla u)-\Delta u\,w\right)dx
=\int_{\partial\Omega}w\,\frac{\partial u}{\partial\nu}\,d\sigma-\int_\Omega\Delta u\,w\,dx
=-\int_\Omega\Delta u\,w\,dx .
$$

（用的是 $w|_{\partial\Omega}=0$。）

**再用带符号的试验函数把 $\Delta u$ 逼为零。** 若存在 $x_0\in\Omega$ 使 $\Delta u(x_0)\ne0$，不妨设 $\Delta u(x_0)>0$。由 $\Delta u$ 的连续性，可取 $\varepsilon>0$ 使 $B_\varepsilon(x_0)\subset\Omega$ 且 $\Delta u$ 在 $B_\varepsilon(x_0)$ 中为正。取

$$
w(x)=\begin{cases}
\bigl(\varepsilon^2-|x-x_0|^2\bigr)^3, & |x-x_0|<\varepsilon,\\[4pt]
0, & |x-x_0|\ge\varepsilon.
\end{cases}
$$

则 $w\in U_0$（在球面上 $w$ 及其一、二阶导数都连续地趋于 $0$），且

$$
\int_\Omega\Delta u\,w\,dx>0,
$$

与 $-\displaystyle\int_\Omega\Delta u\,w\,dx=0$ 矛盾。因此 $\Delta u=0$。

**（$\Leftarrow$）设 $u\in U$ 且 $\Delta u=0$，证明 $u$ 是极小值点。** 任取 $w\in U_0$。在内缩区域上分部积分后取极限，由 $w|_{\partial\Omega}=0$ 得

$$
\left.\frac{d}{dt}E(u+tw)\right|_{t=0}=\int_\Omega\langle\nabla u,\nabla w\rangle\,dx
=\int_{\partial\Omega}w\frac{\partial u}{\partial\nu}\,d\sigma-\int_\Omega\Delta u\,w\,dx=0,
$$

又有

$$
\frac{d^2}{dt^2}E(u+tw)=\int_\Omega\sum_{i=1}^n w_{x_i}^2\,dx\ge0.
$$

因此 $E(u+tw)$ 作为 $t$ 的函数在 $t=0$ 处取极小值。取 $t=1$，即得 $E(u+w)\ge E(u)$，所以 $u$ 为 $E$ 在 $U$ 中的极小值点。$\square$

> **补充（更干净的一步）**　上面的结论可以直接算出来。把被积函数展开：

$$
|\nabla(u+tw)|^2=|\nabla u|^2+2t\langle\nabla u,\nabla w\rangle+t^2|\nabla w|^2,
$$

于是

$$
E(u+tw)=E(u)+t\int_\Omega\langle\nabla u,\nabla w\rangle\,dx+t^2E(w).
$$

只要 $\Delta u=0$（配合 $w|_{\partial\Omega}=0$），一次分部积分给出中间项为 $0$，故

$$
E(u+tw)=E(u)+t^2E(w)\ge E(u)\qquad\text{对一切 }t\in\mathbb R .
$$

这同时说明：$E$ 在 $U$ 上是**凸**泛函，且只要 $w\not\equiv0$ 就有 $E(w)>0$，即严格凸——因此极小元若存在必唯一。

> **补充（Dirichlet 原理为什么重要）**　它把「解 Laplace 方程」换成了「在 $U$ 上求 $E$ 的最小值」。好处是：$E$ 只含**一阶**导数，比 $\Delta u$ 少一阶，于是求解可以在更弱的函数空间里进行；坏处是必须先在某个函数类中证明最小值**取到**（下半连续性 + 强制性），这正是 20 世纪变分法的入口。另一种边界条件下泛函要加边界项，例如

$$
J(u)=\int_\Omega\left(\frac12|\nabla u|^2-fu\right)dx+\int_{\partial\Omega}\left(\frac{\sigma}{2}u^2-gu\right)dS
$$

对应 Robin 边值问题，见 [[12-习题课]] 第一次作业第 4 题。

> **补充（关于 $U$ 与 $U_0$）**　$a+U_0=U$ 是变分法里最常用的一句：它把「在非线性的 $U$ 上做变分」化为「在**线性空间** $U_0$ 上做变分」。后面 [[06-Bochner方法与内部梯度估计]] 里对截断函数做扰动，用的也是同一模式。

---

## 订正

1. **定理 2.1 的陈述有冗余。** 讲稿写作「$u\in U$ 为 $E$ 在 $U$ 中的极小值点，当且仅当 $-\Delta u=0$ 在 $\Omega$ 中、$u=g$ 在 $\partial\Omega$ 上」。但 $U$ 的定义里已经含 $u|_{\partial\Omega}=g$，第二个条件自动成立。等价的说法是：**$u\in U$ 为极小值点 $\iff$ 在 $\Omega$ 内 $\Delta u=0$**。保持原陈述不影响对错，只是稍显重复。
2. **（$\Leftarrow$）方向的关键一步未点明。** 讲稿由「$\left.\frac{d}{dt}E\right|_{t=0}=0$ 且 $\left.\frac{d^2}{dt^2}E\right|_{t=0}\ge0$」直接断言「$t=0$ 处取极小值」。仅凭一阶导为零、二阶导非负是**不够**的（例如 $f(t)=t^3$ 在 $0$ 处同样满足这两条却不是极小值点）。这里之所以成立，是因为 $E(u+tw)$ 关于 $t$ 恰好是**二次多项式**——这一点讲稿没有点出。本篇已在**补充**中用 $E(u+tw)=E(u)+t^2E(w)$ 的形式补全。
3. §2.3.3 的反例中 $\nu(x)=x$ 是因为在单位球面上外法向量就是位置向量，讲稿未标出这一步；本篇已补。

---

系列：[[01-从万有引力到Poisson方程]] · [[03-调和函数·平均值性质与极值原理]] · [[04-调和函数·Harnack不等式与内部梯度估计]]
