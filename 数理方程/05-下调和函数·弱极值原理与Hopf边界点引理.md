# 05 下调和函数·弱极值原理与 Hopf 边界点引理

> 系列：[[04-调和函数·Harnack不等式与内部梯度估计]] · [[06-Bochner方法与内部梯度估计]] · [[07-基本解与广义函数]]

从本篇起进入「解方程」阶段。出发点是把调和函数的不等式方向放宽成**单向**：不再要求 $\Delta u=0$，只要求 $\Delta u\ge0$ 或 $\Delta u\le0$。这样的函数更容易构造，而它们的极值原理却足够强，能把 [[02-边值问题的开始·三类边界条件与Dirichlet原理]] 中留下的 Neumann 问题唯一性补完。

---

## 5.1 定义

**定义**　称 $u\in C^2(\Omega)$ 为**下调和函数**（subharmonic function），若

$$
\Delta u(x)\ge0 .
$$

称 $u\in C^2(\Omega)$ 为**上调和函数**（superharmonic function），若

$$
\Delta u(x)\le0 .
$$

下调和函数很容易实现：

$$
\Delta|x|^2=2n>0 .
$$

**直观**　下调和函数是「某种比凸性弱、但又类似于凸性的东西」：凸函数满足 $f''\ge0$，下调和函数则满足「$\Delta\ge0$」，而 $\Delta$ 是 Hessian 的迹——比逐方向凸要弱。

> **补充（下调和函数 ⟺ 平均值不等式）**　$u$ 下调和当且仅当对一切 $B_r(x_0)\Subset\Omega$ 有

$$
u(x_0)\le\frac{1}{|\partial B_r(x_0)|}\int_{\partial B_r(x_0)}u(x)\,dS_x .
$$

证明完全平行于 [[03-调和函数·平均值性质与极值原理]] 的平均值性质（区别只在于走完最后一步时用 $\Delta u\ge0$ 取代 $\Delta u=0$），完整推导见 [[12-习题课]] 第二次作业第 4 题。因此「下调和」也可以定义成「满足单向平均值性质」。

---

## 5.2 弱极值原理（极大值原理）

**定理（弱极值原理）**　设 $u\in C^2(B_1)\cap C(\overline{B_1})$ 满足 $\Delta u\ge0$，则

$$
\max_{\overline{B_1}}u=\max_{\partial B_1}u .
$$

**证明**　关键是下面这个线性代数事实：**若 $f\in C^2(\Omega)$ 在内点 $x_0\in\Omega$ 处取到极大值，则**

$$
\operatorname{Hess}(f)(x_0)\le0 .
$$

（因为在极大值处沿每个方向的二阶导都非正，故 Hessian 负半定。）

任取 $\varepsilon>0$，令

$$
u_\varepsilon(x):=u+\varepsilon x_1^2 ,
$$

则

$$
\Delta u_\varepsilon=\Delta u+2\varepsilon>0 .
$$

假设 $u_\varepsilon$ 在某个内点 $x_0\in B_1$ 处取到最大值，则

$$
\nabla^2u_\varepsilon(x_0)\le0\quad\Longrightarrow\quad\Delta u_\varepsilon(x_0)\le0,
$$

（$\Delta$ 是 Hessian 的迹，负半定矩阵的迹非正），与 $\Delta u_\varepsilon>0$ 矛盾。**因此 $u_\varepsilon$ 的最大值只能在 $\partial B_1$ 上取到。** 于是对任意 $x\in\overline{B_1}$，

$$
u(x)\le u(x)+\varepsilon x_1^2=u_\varepsilon(x)\le\max_{\partial B_1}u_\varepsilon=\max_{\partial B_1}\bigl(u+\varepsilon x_1^2\bigr)\le\max_{\partial B_1}u+\varepsilon .
$$

由于 $\varepsilon>0$ 任取，令 $\varepsilon\to0^+$ 得

$$
\max_{\overline{B_1}}u\le\max_{\partial B_1}u .
$$

反方向的不等式由 $\partial B_1\subset\overline{B_1}$ 平凡成立。$\square$

**上调和函数的最小值原理**　若 $\Delta u\le0$，对 $-u$ 用上述定理（$\Delta(-u)\ge0$）即得

$$
\min_{\overline{B_1}}u=\min_{\partial B_1}u .
$$

> **补充（对一般有界区域也成立）**　讲稿只对 $B_1$ 陈述，但证明可以逐字搬过去：把扰动函数换成 $u_\varepsilon=u+\varepsilon|x|^2$，则 $\Delta u_\varepsilon=\Delta u+2n\varepsilon>0$，其余论证不变。故对任意有界区域 $\Omega$ 都有 $\max_{\overline\Omega}u=\max_{\partial\Omega}u$。**这一点在下面证明 Hopf 引理时要用于环形区域 $B_1\setminus\overline{B_{3/4}}$。**

> **要点**　弱极值原理不允许 $\varepsilon=0$，因为只假设 $\Delta u\ge0$（不严格）。技巧是「先把方程扰动成严格的、取极限把扰动去掉」。这个手法在后面 [[06-Bochner方法与内部梯度估计]] 里还会再用一次（那里的扰动是「加一项 $C_n u^2$」）。

---

## 5.3 强极值原理（Hopf 边界点引理）

弱极值原理只说明「极大值在边界取到」，但**没有说极大值点的法向导数会怎样**。强极值原理补上了这一点。

**定理（Hopf 边界点引理）**　设 $u\in C^2(B_1)\cap C(\overline{B_1})$ 满足 $\Delta u\ge0$。若 $x_0\in\partial B_1$ 处满足

$$
u(x)<u(x_0),\qquad\forall x\in B_1,
$$

则

$$
\liminf_{t\to0^+}\frac{u(x_0)-u(x_0-t\nu)}{t}>0 .
$$

其中 $\nu$ 为 $x_0$ 处单位外法向量。若 $u$ 在 $x_0$ 处可微，则

$$
\frac{\partial u}{\partial\nu}(x_0)>0 .
$$

**三条观察**

1. **例**：$|x|^2$ 在 $\partial B_R$ 上满足

$$
\frac{\partial|x|^2}{\partial\nu}\Big|_{\partial B_R}
=\bigl\langle 2x,\tfrac{x}{|x|}\bigr\rangle\Big|_{|x|=R}=2R>0 .
$$

2. 在上述引理的假设下，若 $u$ 在 $x_0$ 处可微，则对**任意与 $\nu$ 成锐角的单位向量** $e$（即 $\langle e,\nu\rangle>0$），有

$$
\frac{\partial u}{\partial e}(x_0)=\langle e,\nu\rangle\frac{\partial u}{\partial\nu}(x_0)>0 .
$$

也就是说：只要「朝外」走，方向导数就严格为正。

3. 若 $u\in C^2(B_1)\cap C(\overline{B_1})$ 满足 $\Delta u\le0$，且在 $x_0\in\partial B_1$ 处满足 $u(x)>u(x_0)$（$\forall x\in B_1$），则对 $-u$ 应用上述引理可知：若 $u$ 在 $x_0$ 处可微，则

$$
\frac{\partial u}{\partial\nu}(x_0)<0 .
$$

### 5.3.1 证明

证明的思想是：把一个已知的、在边界上法向导数为负的**障碍函数**（barrier）加到 $u$ 上去，再用**弱**极值原理。

**障碍函数**　记 $\Gamma_n$ 为基本解（见 [[07-基本解与广义函数]]）：

$$
\Gamma_n(x)=\begin{cases}
-\dfrac{1}{2\pi}\log|x|, & n=2,\\[8pt]
\dfrac{1}{(n-2)|S^{n-1}|\,|x|^{n-2}}, & n\ge3 .
\end{cases}
$$

它有三条性质：$\Gamma_n\ge0$ 于 $B_1\setminus\{0\}$，$\Delta\Gamma_n=0$ 于 $\mathbb R^n\setminus\{0\}$，且

$$
\frac{\partial\Gamma_n}{\partial\nu}\Big|_{\partial B_1}<0 .
$$

> **补充（第三条性质的验算）**　$n\ge3$ 时 $\nabla\Gamma_n(x)=-\dfrac{x}{|S^{n-1}||x|^{n}}$，故在 $\partial B_1$ 上

$$
\frac{\partial\Gamma_n}{\partial\nu}=\bigl\langle\nabla\Gamma_n,x\bigr\rangle=-\frac{1}{|S^{n-1}|}<0 ;
$$

> $n=2$ 时 $\nabla\Gamma_2(x)=-\dfrac{1}{2\pi}\dfrac{x}{|x|^2}$，故在 $\partial B_1$ 上 $\dfrac{\partial\Gamma_2}{\partial\nu}=-\dfrac{1}{2\pi}<0$。同时 $\Gamma_2(x)=-\dfrac{1}{2\pi}\log|x|>0$ 当 $0<|x|<1$。

**证明**　设 $u\in C^2(B_1)\cap C(\overline{B_1})$ 满足 $\Delta u\ge0$，且 $u$ 在 $x_0\in\partial B_1$ 处满足

$$
u(x_0)>u(x),\qquad x\in B_1,
$$

则（由紧性与连续性）

$$
u(x_0)>\max_{\partial B_{3/4}}u .
$$

取 $\varepsilon>0$ 充分小，使得

$$
\max_{\partial B_{3/4}}\bigl(u(x)+\varepsilon\Gamma_n(x)\bigr)<u(x_0).
$$

（这能做到，因为 $\Gamma_n$ 在 $\partial B_{3/4}$ 上有界，而 $\varepsilon$ 可以任意小。）

**在环域上用弱极值原理。** 记 $v=u+\varepsilon\Gamma_n$。则

$$
\Delta v=\Delta u+\varepsilon\Delta\Gamma_n=\Delta u\ge0\qquad\text{于 }B_1\setminus\overline{B_{3/4}},
$$

即 $v$ 在环域上下调和。又 $\Gamma_n$ 在 $\partial B_1$ 上等于非负常数 $\Gamma_n(x_0)$，而 $u\le u(x_0)$ 于 $\overline{B_1}$，故 $v$ 在环域的外边界 $\partial B_1$ 上满足

$$
v\le u(x_0)+\varepsilon\Gamma_n(x_0)=v(x_0).
$$

在内边界 $\partial B_{3/4}$ 上由 $\varepsilon$ 的取法有 $v<u(x_0)\le v(x_0)$。于是在环域的**两条**边界上都有 $v\le v(x_0)$。由弱极值原理（对一般有界区域的形式，见 5.2 的补充），环域内处处有

$$
v\le v(x_0).
$$

因此对 $t>0$，

$$
v(x_0)-v(x_0-t\nu)\ge0,
$$

即

$$
\liminf_{t\to0^+}\frac{v(x_0)-v(x_0-t\nu)}{t}\ge0 .
$$

**把 $\varepsilon\Gamma_n$ 的部分剥掉。** 由 $v=u+\varepsilon\Gamma_n$，

$$
\liminf_{t\to0^+}\frac{u(x_0)-u(x_0-t\nu)}{t}
=\liminf_{t\to0^+}\frac{v(x_0)-v(x_0-t\nu)}{t}
-\lim_{t\to0^+}\frac{\varepsilon\Gamma_n(x_0)-\varepsilon\Gamma_n(x_0-t\nu)}{t}.
$$

第二个极限是 $\varepsilon$ 乘 $\Gamma_n$ 在 $x_0$ 处沿 $\nu$ 的方向导数，即 $\varepsilon\dfrac{\partial\Gamma_n}{\partial\nu}(x_0)$。于是

$$
\liminf_{t\to0^+}\frac{u(x_0)-u(x_0-t\nu)}{t}\ge-\varepsilon\frac{\partial\Gamma_n}{\partial\nu}(x_0)>0 ,
$$

最后一步用了 $\dfrac{\partial\Gamma_n}{\partial\nu}(x_0)<0$。$\square$

**推广到一般区域**　对于具有 $C^2$ 边界的区域，可在边界点处取**内切球**，再经平移、缩放，把上述单位球上的 Hopf 引理用于该球。内切球与区域在切点处的外法向量相同，故结论照搬。

> **要点**　证明的机制值得记牢：**障碍函数 $\Gamma_n$ 的作用是「提供一个在边界上法向严格为负的已知量」**。因为要证明的量是 $\liminf$（不是函数值本身），所以可以把它拆成「$v$ 的差商」减去「$\varepsilon\Gamma_n$ 的差商」；前者由弱极值原理非负，后者由 $\partial\Gamma_n/\partial\nu<0$ 严格为负，相减得到严格为正。这里 $u\ge0$ 完全没用到——**Hopf 引理不需要 $u$ 非负**，比 [[04-调和函数·Harnack不等式与内部梯度估计]] 的梯度估计假设更弱。

---

## 5.4 两个推论：Neumann 问题解的唯一性

**推论 1（Neumann 内问题）**　设 $n\ge2$，$\Omega\subset\mathbb R^n$ 为有界连通域，边界为 $C^2$。对固定的 $f$ 和 $g$，若

$$
\begin{cases}
\Delta u=f, & \text{于 }\Omega,\\[4pt]
\dfrac{\partial u}{\partial\nu}=g, & \text{于 }\partial\Omega
\end{cases}
$$

在 $C^2(\Omega)\cap C^1(\overline\Omega)$ 中有解，则**任意两个解相差一个常数**。这里 $\nu$ 为 $\Omega$ 的单位外法向量。

**证明**　设 $u_1,u_2$ 为两个解，令 $w=u_1-u_2$，则

$$
\Delta w=0\ \text{于 }\Omega,\qquad \frac{\partial w}{\partial\nu}=0\ \text{于 }\partial\Omega .
$$

假设 $w$ 不为常数。由弱极值原理，$w$ 的最大值只能在边界取到，故存在 $x_0\in\partial\Omega$ 使得

$$
w(x)<w(x_0),\qquad x\in\Omega
$$

（严格不等号来自「$w$ 非常数 + 强极值原理」，即最大值不在内部、也不能在内部取到与边界相同的值）。由 Hopf 引理（用于 $x_0$ 处的内切球），

$$
\frac{\partial w}{\partial\nu}(x_0)>0,
$$

与齐次边界条件 $\dfrac{\partial w}{\partial\nu}=0$ 矛盾。因此 $w$ 为常数。$\square$

**推论 2（Neumann 外问题）**　设 $n\ge2$，$\Omega\subset\mathbb R^n$ 为具有 $C^2$ 边界的有界域，且 $E:=\mathbb R^n\setminus\overline\Omega$ 连通。对固定的 $f$ 和 $g$，若

$$
\begin{cases}
\Delta u=f, & \text{于 }E,\\[4pt]
\dfrac{\partial u}{\partial\nu}=g, & \text{于 }\partial E,\\[4pt]
\lim\limits_{|x|\to\infty}u(x)=0
\end{cases}
$$

在 $C^2(E)\cap C^1_{\mathrm{loc}}(E)$ 中有解，则**解唯一**。这里 $\nu$ 为 $E$ 的单位外法向量，即在 $\partial E$ 上指向 $\Omega$。

**证明**　设 $u_1,u_2$ 为两个解，令 $w=u_1-u_2$，则

$$
\Delta w=0\ \text{于 }E,\qquad \frac{\partial w}{\partial\nu}=0\ \text{于 }\partial E,\qquad \lim_{|x|\to\infty}w(x)=0 .
$$

假设 $w\not\equiv0$。必要时将 $w$ 换为 $-w$，可设存在 $y\in E$ 使得 $w(y)>0$。取 $R>|y|$ 充分大，使得

$$
|w(x)|<\frac{w(y)}{2},\qquad x\in E\setminus B_R .
$$

于是 $w$ 在 $E$ 上的**正最大值**在某点 $x_0$ 处取得（因为无穷远处 $w$ 已经小于 $w(y)/2$，而 $y$ 处取到正値，故正最大值不会跑到无穷远）。由强极值原理，$x_0$ 不可能位于 $E$ 的内部（否则 $w$ 在 $E$ 上恒为常数，与 $\lim w=0$ 且 $w(y)>0$ 矛盾），故 $x_0\in\partial E$，且

$$
w(x)<w(x_0),\qquad x\in E .
$$

由 Hopf 引理，$\dfrac{\partial w}{\partial\nu}(x_0)>0$，与齐次 Neumann 条件矛盾。因此 $w\equiv0$，即 $u_1=u_2$。$\square$

> **补充（为什么必须用强极值原理）**　只用弱极值原理只能得到「$w$ 的最大值在边界取到」，但这与 $\partial w/\partial\nu=0$ 并不矛盾——最大值在边界取到完全可以是「平的」。Hopf 引理的价值恰恰在于：**只要非常数，边界极大值点处的法向导数就严格为正**，于是与齐次 Neumann 条件直接冲突。这就是强极值原理在 Neumann 问题里不可替代的原因。

> **补充（Neumann 内问题何时真有解）**　推论 1 只谈「若有解则相差常数」。有解本身需要相容性条件：对 $\Delta u=f$ 在 $\Omega$ 上积分并用散度定理得

$$
\int_\Omega f\,dx=\int_{\partial\Omega}\frac{\partial u}{\partial\nu}\,d\sigma=\int_{\partial\Omega}g\,d\sigma ,
$$

即 $\displaystyle\int_\Omega f=\int_{\partial\Omega}g$ 是必要条件（对任一解成立）。这与「相差常数」并不冲突：常数可以自由调整，但相容性条件与常数无关。

---

## 订正

1. **基本解在二维的符号，前后两讲不一致。** 讲稿第二周开头把二维的核写成

$$
G(x,y)=\log|x-y|\quad(n=2),\qquad G(x,y)=\frac{1}{|x-y|^{n-2}}\quad(n\ge3);
$$

而第三周起统一写成

$$
\Gamma_2(x)=-\frac{1}{2\pi}\log|x|,\qquad \Gamma_n(x)=\frac{1}{(n-2)|S^{n-1}||x|^{n-2}}\ (n\ge3).
$$

两者在 $n=2$ 时**相差一个整体符号**。由于第二周那处只用来验证「$x\ne y$ 时调和」（$\pm\log|x|$ 都调和），结论不受影响；但作为**基本解**的约定必须统一为带负号的那个，否则 $-\Delta\Gamma=\delta_0$ 与 $\partial\Gamma/\partial\nu\big|_{\partial B_1}=-1<0$ 都会反号。本篇与 [[07-基本解与广义函数]]、[[08-Green恒等式与Green表示公式]]、[[09-Green函数的性质]] 一律采用 $-\dfrac{1}{2\pi}\log|x|$。
2. **弱极值原理的光滑性假设，讲稿两处写法不同。** 第二周写 $u\in C^2(B_1)\cap C^1(\overline{B_1})$，第三周写 $u\in C^2(B_1)\cap C(\overline{B_1})$。证明中只用到 $u$ 在 $\overline{B_1}$ 上连续（保证最大值取到）与在 $B_1$ 内 $C^2$（保证 Hessian 存在），故 $C(\overline{B_1})$ 是正确而更弱的假设，本篇采用后者。
3. **Hopf 引理中极限的写法，讲稿两处不同。** 第二周写 $\lim_{t\to0^+}$，第三周写 $\liminf_{t\to0^+}$。证明中用到的是「$v(x_0)-v(x_0-t\nu)\ge0$ 对一切充分小 $t>0$」，因此结论天然是 $\liminf$；在 $u\in C^2(\overline{B_1})$（极限存在）时两者等价。本篇统一写作 $\liminf$（更弱、更保险）。
4. **「$w$ 的最大值只能在边界取到」到「$w(x)<w(x_0)$ 对一切 $x\in\Omega$」这一步讲稿未说明。** 这一步恰恰是**强**极值原理的作用（非常数的上调和/调和函数在内部不能取到最大值），不能只靠弱极值原理推出。本篇已在推论 1 的证明中补上括号说明。
5. **第三周 §2 中「$u+\varepsilon\Gamma_n$ 在环域上下调和」的直接原因是 $\Delta\Gamma_n=0$**（而不只是「$\Gamma_n$ 调和」这一说法本身）。本篇已写明。

---

系列：[[04-调和函数·Harnack不等式与内部梯度估计]] · [[06-Bochner方法与内部梯度估计]] · [[07-基本解与广义函数]]
