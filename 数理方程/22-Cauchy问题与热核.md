# 22 Cauchy 问题与热核

> 系列：[[21-Fourier变换与速降函数空间]] · [[23-第三次作业]] · [[19-分离变量解的正则性·三个收敛命题]]

前面三篇都在有界区域 $\Omega$ 里做事：分离变量、特征展开、Duhamel 叠加。这一篇把 $\Omega$ 换成整个 $\mathbb R^n$，问同一个热方程还怎么解。

答案只有一句：**把 [[21-Fourier变换与速降函数空间]] 的 Fourier 变换套上去，PDE 变成 ODE，解出来再变回来。** 变回来的那个核就是**热核** $G_t(x)=\dfrac{1}{(4\pi t)^{n/2}}e^{-\frac{|x|^2}{4t}}$。

讲稿名义上叫「Cauchy 问题与热核」，实际上热核的积分只算到配方那一步（第 12 页末）就停了 —— 常数因子、换元、以及「这东西为什么真的是热核」全没写。这一篇把它补完，并把讲稿那个双重计数的排印错误指出来（订正第 2 条）。

---

## 22.1 问题的提出：全空间里「分离变量」失效

**Cauchy 问题**（全空间热方程，初值问题）

$$
\begin{cases}\dfrac{\partial u}{\partial t}=\Delta u, & (x,t)\in\mathbb R^n\times(0,+\infty),\\[6pt] u(x,0)=\varphi(x), & x\in\mathbb R^n .\end{cases}\tag{22.1}
$$

与 [[18-热方程·建模、边值条件与分离变量法]] 的初边值问题相比，**只少了边界条件**。但正是这一「少」把整个方法换掉了：

| | 有界区域 $\Omega$ | 全空间 $\mathbb R^n$ |
|---|---|---|
| 谱 | 离散 $\lambda_k\to+\infty$ | 连续谱 $[0,+\infty)$（§[[21-Fourier变换与速降函数空间]] 21.1） |
| 展开 | 级数 $\sum_k a_ke^{-\lambda_kt}u_k(x)$ | 积分（逆 Fourier 变换） |
| 边界的角色 | 施加 Dirichlet/Neumann/Robin，制造离散谱 | **没有边界**，改要求 $\varphi\in\mathcal S(\mathbb R^n)$ 保证积分收敛 |
| 解的核 | $p_\Omega(x,y;t)=\sum_ke^{-\lambda_kt}u_k(x)u_k(y)$ | $G_t(x-y)$ |

最后一行是两边的对应点：有界区域的解 $u(x,t)=\int_\Omega p_\Omega(x,y;t)\varphi(y)dy$，全空间的解 $u(x,t)=\int_{\mathbb R^n}G_t(x-y)\varphi(y)dy$。**「热核」就是全空间的那个 $p$。**

> **要点**　$p_\Omega$ 与 $G_t$ 的关系后面还会用到：$p_\Omega$ 是 $G_t$ 在 $\Omega$ 上「用特征函数做的谱投影」，而 $G_t$ 可以看成 $\Omega=\mathbb R^n$ 时 $p$ 的极限（此时 $\sum_k$ 变成 $\int d\xi$，特征函数变成 $e^{\sqrt{-1}\langle\xi,x\rangle}$）。见 §22.8。

---

## 22.2 Fourier 侧：PDE 变成 ODE

对 (22.1) 两边关于 $x$ 作 Fourier 变换。由 [[21-Fourier变换与速降函数空间]] 命题 21.4.1，

$$
\mathcal F[\Delta u(\cdot,t)](\xi)=\sum_{j=1}^{n}(\sqrt{-1}\xi_j)^2\widehat u(\xi,t)=-|\xi|^2\widehat u(\xi,t),
$$

而关于 $t$ 的导数与 $\mathcal F$ 可交换（$\mathcal F$ 只作用在 $x$ 上；且 $u$ 速降使积分号下求导合法，见订正第 4 条）：

$$
\mathcal F\left[\frac{\partial u}{\partial t}(\cdot,t)\right](\xi)=\frac{d}{dt}\widehat u(\xi,t).
$$

于是 (22.1) 变成一族**常微分方程**（每个 $\xi$ 一个）：

$$
\begin{cases}\dfrac{d}{dt}\widehat u(\xi,t)=-|\xi|^2\widehat u(\xi,t), & t>0,\\[8pt] \widehat u(\xi,0)=\widehat\varphi(\xi).\end{cases}\tag{22.2}
$$

这是最平凡的一阶线性方程，解为

$$
\boxed{\ \widehat u(\xi,t)=e^{-|\xi|^2t}\,\widehat\varphi(\xi)\ }\tag{22.3}
$$

> **要点（$e^{-|\xi|^2t}$ 就是「衰减因子」）**　$(\sqrt{-1}\xi_j)^2=-\xi_j^2$ 这个负号是热方程区别于波方程的全部原因：它让 (22.2) 的解是指数**衰减**而不是振荡。$|\xi|$ 越大衰减越快 —— **高频分量在任意小的正时刻就被抹掉了**，这就是 §22.6 里「$t>0$ 立刻无穷次可微」的机制。对照 [[19-分离变量解的正则性·三个收敛命题]]，那里的 $e^{-\lambda_kt}$ 与这里的 $e^{-|\xi|^2t}$ 是同一件事，只是 $\lambda_k$ 离散、$|\xi|^2$ 连续。

---

## 22.3 从 Fourier 侧回到物理侧：卷积

(22.3) 右端是**两个函数的乘积**。由 [[21-Fourier变换与速降函数空间]] 定理 21.4.2 的卷积定理，乘积的逆变换是卷积的逆变换：

$$
u(\cdot,t)=\mathcal F^{-1}\left[e^{-|\xi|^2t}\widehat\varphi\right]=G_t*\varphi ,\qquad G_t:=\mathcal F^{-1}\left[e^{-|\xi|^2t}\right].\tag{22.4}
$$

即

$$
u(x,t)=\int_{\mathbb R^n}G_t(x-y)\varphi(y)\,dy .\tag{22.5}
$$

所以**整个问题归结为一个积分**：算 $G_t$。这个 $G_t$ 就叫 $n$ 维**热核**（heat kernel），也叫 Gauss–Weierstrass 核。

---

## 22.4 算热核

由定义 21.3.1 的逆变换公式（注意那个 $\frac{1}{(2\pi)^n}$，讲稿漏了，见订正第 1 条），

$$
G_t(x)=\mathcal F^{-1}\left[e^{-|\xi|^2t}\right](x)=\frac{1}{(2\pi)^n}\int_{\mathbb R^n}e^{-t|\xi|^2+\sqrt{-1}\langle\xi,x\rangle}\,d\xi .
$$

### 22.4.1 配方

把指数里对 $\xi$ 的二次式配成完全平方。目标是写成 $|\xi-c|^2$ 的形状，$c$ 待定：

$$
|\xi-c|^2=|\xi|^2-2\langle\xi,c\rangle+|c|^2 .
$$

要让 $-2\langle\xi,c\rangle$ 对上 $+\dfrac{\sqrt{-1}}{t}\langle\xi,x\rangle$，取

$$
c=\frac{\sqrt{-1}\,x}{2t}\qquad\Longrightarrow\qquad |c|^2=\frac{(\sqrt{-1})^2|x|^2}{4t^2}=-\frac{|x|^2}{4t^2}.
$$

于是

$$
|\xi|^2-\frac{\sqrt{-1}\langle\xi,x\rangle}{t}=\left|\xi-\frac{\sqrt{-1}x}{2t}\right|^2+\frac{|x|^2}{4t^2},
$$

两边乘 $-t$：

$$
-t|\xi|^2+\sqrt{-1}\langle\xi,x\rangle=-t\left|\xi-\frac{\sqrt{-1}x}{2t}\right|^2-\frac{|x|^2}{4t} .\tag{22.6}
$$

**这就是那个 $e^{-\frac{|x|^2}{4t}}$ 的来源 —— 它来自配方里的 $|c|^2$，只能出现一次。** 讲稿把它在括号内（$+\frac{|x|^2}{4t^2}$ 乘 $-t$ 得 $-\frac{|x|^2}{4t}$）和括号外各算了一遍，见订正第 2 条。

### 22.4.2 复数平移

把 (22.6) 代回：

$$
G_t(x)=\frac{1}{(2\pi)^n}\left[\int_{\mathbb R^n}e^{-t\left|\xi-\frac{\sqrt{-1}x}{2t}\right|^2}d\xi\right]e^{-\frac{|x|^2}{4t}} .
$$

接下来换元 $\eta=\xi-\dfrac{\sqrt{-1}x}{2t}$。**这是全篇唯一不平凡的一步**：平移量是纯虚的，不是实平移，不能直接说「换元不变」。合法性如下（$n=1$ 时就是 Cauchy 定理）：

固定 $x$、$t>0$，令 $F(\zeta):=e^{-t\zeta^2}$（$\zeta\in\mathbb C$）为整函数。对矩形围道 $[-R,R]\to[R,R-\sqrt{-1}x_1/2t]\to[-R,-R-\sqrt{-1}x_1/2t]$，$F$ 在其上解析，故沿闭围道积分为 $0$。两条竖直边的积分在 $R\to\infty$ 时趋于 $0$，因为 $\zeta=\pm R+\sqrt{-1}s$ 上

$$
\left|e^{-t\zeta^2}\right|=e^{-t(R^2-s^2)}=e^{-tR^2}e^{ts^2},
$$

而 $s$ 只在有界区间 $[0,-x_1/2t]$ 里跑，故被 $e^{-tR^2+ts_{\max}^2}\to0$ 压住。于是水平边上的积分相等，即

$$
\int_{\mathbb R}e^{-t\xi_1^2}d\xi_1=\int_{\mathbb R}e^{-t\left(\xi_1-\frac{\sqrt{-1}x_1}{2t}\right)^2}d\xi_1 .
$$

$n$ 维逐变量做一遍（被积函数是各变量的一元整函数之积，其余变量固定），结果一样。

### 22.4.3 Gauss 积分

平移之后剩下的是一元 Gauss 积分：

$$
\int_{\mathbb R}e^{-t\eta_j^2}d\eta_j=\sqrt{\frac{\pi}{t}},
$$

故

$$
\int_{\mathbb R^n}e^{-t|\eta|^2}d\eta=\left(\frac{\pi}{t}\right)^{n/2}.
$$

### 22.4.4 结果

$$
\boxed{\ G_t(x)=\frac{1}{(2\pi)^n}\left(\frac{\pi}{t}\right)^{n/2}e^{-\frac{|x|^2}{4t}}=\frac{1}{(4\pi t)^{n/2}}\,e^{-\frac{|x|^2}{4t}}\ }\tag{22.7}
$$

（因为 $(4\pi t)^{n/2}=2^n\pi^{n/2}t^{n/2}=(2\pi)^n(\pi/t)^{n/2}$。）

一维、三维写开：

$$
G_t(x)=\frac{1}{\sqrt{4\pi t}}e^{-\frac{x^2}{4t}}\ (n=1),\qquad G_t(x)=\frac{1}{(4\pi t)^{3/2}}e^{-\frac{|x|^2}{4t}}\ (n=3).
$$

> **补充（$(4\pi t)^{n/2}$ 是被 $\int G_t=1$ 逼出来的）**　这个常数没有任何自由：任何满足 (22.1) 且「质量守恒」的解都要求 $G_t$ 积分为 $1$（§22.5 第 1 条）。若照讲稿那个双重计数的式子算，会得到 $\frac{1}{(2\pi)^n}(\pi/t)^{n/2}e^{-\frac{|x|^2}{2t}}$，积分为 $\frac{1}{(2\pi)^n}(\pi/t)^{n/2}(2\pi t)^{n/2}=2^{-n}\ne1$ —— **这就是检测那个排印错误的最快办法。** 另外，若改用 [[21-Fourier变换与速降函数空间]] §21.3 表里的「对称约定」，则 $\mathcal F[G_t]=(2\pi)^{-n/2}\!\int G_te^{-\sqrt{-1}\langle x,\xi\rangle}dx$，同样得到 (22.7)；换约定时 $2\pi$ 的账在正变换与逆变换之间来回搬，核本身不变。

---

## 22.5 热核的基本性质

**1. 正性与单位质量**

$$
G_t(x)>0\quad(x\in\mathbb R^n),\qquad \int_{\mathbb R^n}G_t(x)\,dx=1\quad(t>0).
$$

后者由 Gauss 积分：$\int_{\mathbb R}e^{-\frac{x_j^2}{4t}}dx_j=\sqrt{4\pi t}$，$n$ 个相乘得 $(4\pi t)^{n/2}$，正好抵掉前面的常数。

质量守恒的物理含义就是「热量不增不减，只是重新分布」。数学上它保证 (22.5) 是 $\varphi$ 的一个**加权平均**：

$$
\min\varphi\le u(x,t)\le\max\varphi ,
$$

即**极值原理在全空间也成立**（对照 [[05-下调和函数·弱极值原理与Hopf边界点引理]]）。

**2. 半群性质**

$$
G_t*G_s=G_{t+s}\qquad(t,s>0).
$$

**证明**　由卷积定理与 (22.7)，

$$
\mathcal F[G_t*G_s](\xi)=e^{-|\xi|^2t}\cdot e^{-|\xi|^2s}=e^{-|\xi|^2(t+s)}=\mathcal F[G_{t+s}](\xi),
$$

而 $G_t*G_s$ 与 $G_{t+s}$ 都在 $\mathcal S$ 中、$\mathcal F$ 在 $\mathcal S$ 上单射（[[21-Fourier变换与速降函数空间]] §21.4.4），故相等。$\square$

这正是 [[20-Duhamel原理与齐次化]] §20.5 里 $e^{tA}e^{sA}=e^{(t+s)A}$ 的积分形式：**「先扩散 $s$ 再扩散 $t$」与「一次扩散 $t+s$」等价。**

**3. 自身满足热方程**

$$
\frac{\partial}{\partial t}G_t(x)=\Delta_xG_t(x)\qquad(t>0).
$$

直接验：由 $G_t(x)=(4\pi t)^{-n/2}e^{-|x|^2/4t}$，

$$
\partial_{x_j}G_t=-\frac{x_j}{2t}G_t,\qquad \partial^2_{x_j}G_t=\left(\frac{x_j^2}{4t^2}-\frac{1}{2t}\right)G_t ,
$$

求和得 $\Delta G_t=\left(\dfrac{|x|^2}{4t^2}-\dfrac{n}{2t}\right)G_t$；而对 $t$ 求导（$e^{-\frac{|x|^2}{4t}}$ 里 $t$ 在分母上，导数带正号）得 $\partial_tG_t=\left(-\dfrac{n}{2t}+\dfrac{|x|^2}{4t^2}\right)G_t$。两者相同。

> **要点**　所以 $G_t$ 本身是 (22.1) 取 $\varphi=\delta_0$ 的解（在广义函数意义下）—— 它是**基本解**（[[07-基本解与广义函数]]）。$n=1$ 时这与 $\Gamma(x)=|x|$ 的关系是：$\Gamma$ 解 $-\Delta\Gamma=\delta_0$（椭圆、稳态），$G_t$ 解 $\partial_tG_t-\Delta G_t=0$ 且 $G_t\to\delta_0$（抛物、演化）。

**4. 集中性（近似单位）**

对任意 $\delta>0$，

$$
\int_{|x|>\delta}G_t(x)\,dx\le 2^{n/2}e^{-\frac{\delta^2}{8t}}\xrightarrow[t\to0^+]{}\ 0 .
$$

**证明**　$|x|>\delta$ 时 $\dfrac{|x|^2}{4t}=\dfrac{|x|^2}{8t}+\dfrac{|x|^2}{8t}>\dfrac{\delta^2}{8t}+\dfrac{|x|^2}{8t}$，故

$$
\int_{|x|>\delta}G_t\,dx\le e^{-\frac{\delta^2}{8t}}\cdot\frac{1}{(4\pi t)^{n/2}}\int_{\mathbb R^n}e^{-\frac{|x|^2}{8t}}dx=e^{-\frac{\delta^2}{8t}}\cdot\frac{(8\pi t)^{n/2}}{(4\pi t)^{n/2}}=2^{n/2}e^{-\frac{\delta^2}{8t}} .
$$

$\square$

即 $t\to0^+$ 时 $G_t$ 的**全部质量向原点收缩**。这是 (22.5) 能在 $t\to0^+$ 时还原出 $\varphi$ 的根据。

**5. 概率解释**

$G_t$ 是 $\mathcal N(0,2tI_n)$ 的密度。设 $(B_t)_{t\ge0}$ 为 $n$ 维标准 Brown 运动，则 $\sqrt{2}B_t\sim\mathcal N(0,2tI_n)$，于是 (22.5) 写成

$$
u(x,t)=\mathbb E\left[\varphi\left(x+\sqrt{2}B_t\right)\right].
$$

这就是 Feynman–Kac 公式在无势情形下的形式：**热方程的解 = 初值沿 Brown 运动的平均。** 立刻可以验两个例子：$\varphi\equiv1$ 给 $u\equiv1$；$\varphi(x)=|x|^2$ 给 $u(x,t)=|x|^2+2nt$（因为 $\mathbb E|B_t|^2=nt$）—— 后者可以直接代进 (22.1) 检查（$\partial_tu=2n$，$\Delta u=2n$），也是个好用的 sanity check。

---

## 22.6 初值的取回与正则化

### 22.6.1 $t>0$ 立刻 $C^\infty$

**命题 22.6.1**　设 $\varphi\in L^\infty(\mathbb R^n)$，$u$ 由 (22.5) 给出。则 $u\in C^\infty(\mathbb R^n\times(0,+\infty))$，且对任意多重指标 $\alpha$、$k\ge0$，

$$
\partial_t^k\partial_x^\alpha u(x,t)=\int_{\mathbb R^n}\partial_t^k\partial_x^\alpha G_t(x-y)\,\varphi(y)\,dy .\tag{22.8}
$$

**证明**　固定 $t_0>0$，在 $t\ge t_0$ 上，对任意 $\alpha,k$ 有

$$
\left|\partial_t^k\partial_x^\alpha G_t(x-y)\right|\le C_{\alpha,k,t_0}\qquad(\text{对一切 }x,y,t\ge t_0),
$$

因为 $G_t$ 及其各阶导数都是「多项式 $\times$ 指数衰减」，且 $t\ge t_0$ 时指数衰减率有一致下界（$\frac{|x-y|^2}{4t}\ge\frac{|x-y|^2}{4T}$，$T$ 为所考察的 $t$ 的上界）。又 $\varphi\in L^\infty$，故被积函数被可积函数 $C_{\alpha,k,t_0}\lVert\varphi\rVert_\infty$ 控制，由控制收敛定理可在积分号下求任意次导数。(22.8) 中的积分对 $(x,t)$ 连续（同样理由）。$\square$

> **要点（抛物型正则化）**　**初值可以很坏，解在任意小的正时刻就变光滑。** 这与 [[19-分离变量解的正则性·三个收敛命题]] 的结论方向一致，但强得多：那里要 $C^2$ 相容性才拿到 $C^1_tC^2_x$，这里 $\varphi$ 只要**有界可测**（不必连续、不必可积）。原因在于那里的光滑性来自「初值的相容性 + 系数 $e^{-\lambda_kt}$ 的衰减」，只能拿到初值已有的那么多阶；这里的光滑性来自 $|\xi|^2\to\infty$ 时 $e^{-|\xi|^2t}$ 的**超强衰减**，它对所有频率一视同仁地压平，所以给出无穷阶。热方程是全门课里唯一有这种性质的方程（波方程没有：$\cos(|\xi|t)$ 不衰减）。

### 22.6.2 $t\to0^+$ 的三种收敛

**命题 22.6.2**　设 $\varphi\in C_b(\mathbb R^n)$（有界连续）。则 $u(\cdot,t)\to\varphi$ 在 $\mathbb R^n$ 上**局部一致**；若 $\varphi$ 还一致连续，则**整体一致**。

**证明**　由 $\int G_t=1$，

$$
u(x,t)-\varphi(x)=\int_{\mathbb R^n}G_t(y)\left[\varphi(x-y)-\varphi(x)\right]dy ,
$$

于是对任意 $\delta>0$，

$$
\left|u(x,t)-\varphi(x)\right|\le\sup_{|y|\le\delta}\left|\varphi(x-y)-\varphi(x)\right|+2\lVert\varphi\rVert_\infty\int_{|y|>\delta}G_t\,dy .
$$

第二项由 §22.5 第 4 条 $\le 2^{n/2+1}\lVert\varphi\rVert_\infty e^{-\frac{\delta^2}{8t}}\to0$（$t\to0^+$）。第一项：在紧集 $K$ 上 $\varphi$ 一致连续，取 $\delta$ 使 $|y|\le\delta$ 时 $\sup_{x\in K}|\varphi(x-y)-\varphi(x)|<\varepsilon$。先定 $\delta$，再取 $t$ 使第二项 $<\varepsilon$。$\square$

若 $\varphi$ 只是**有界可测**，则同一估计给出：在 $\varphi$ 的每个**连续点** $x$ 处 $u(x,t)\to\varphi(x)$；若 $\varphi\in L^p$（$1\le p<\infty$），则 $u(\cdot,t)\to\varphi$ 于 $L^p$。特别地 $\varphi\in L^2$ 时由 [[21-Fourier变换与速降函数空间]] 定理 21.5.1 与 (22.3)，

$$
\lVert u(\cdot,t)-\varphi\rVert_{L^2}^2=\frac{1}{(2\pi)^n}\int_{\mathbb R^n}\left|\widehat\varphi(\xi)\right|^2\left|e^{-|\xi|^2t}-1\right|^2d\xi\xrightarrow[t\to0^+]{}0 ,
$$

因为 $|e^{-|\xi|^2t}-1|\le1$ 且逐点趋于 $0$，用控制收敛（被 $|\widehat\varphi|^2$ 控制，它在 $L^1$ 里）。

三种收敛与 [[19-分离变量解的正则性·三个收敛命题]] 的对照：

| | 有界区域（19 篇） | 全空间（本篇） |
|---|---|---|
| 内闭一致 + 各阶导数 | 需 $C^2$ 相容性 | **不需要**任何相容性（没有边界） |
| $L^2$ 收敛 | 需 $\sum|a_k|^2k^2<\infty$ 一类 | 需 $\varphi\in L^2$ |
| $C^0$ 收敛 | 需 $C^0$ 相容性（角点连续） | 需 $\varphi$ 在 $x$ 处连续 |

差别全在「有没有边界」：**有界区域的不光滑性来自边界角点的失配，全空间的不光滑性只来自 $\varphi$ 自己。**

### 22.6.3 唯一性：为什么必须限制增长

**补充（Tychonoff 反例）**　在 $\mathcal S$ 类里，(22.3) 给出的是唯一解（$\mathcal F$ 单射，加 (22.2) 的 ODE 唯一性）。但**一旦不限制增长，(22.1) 的唯一性就失效**：Tychonoff（1935）构造了 $\mathbb R^n\times(0,T)$ 上 $u\not\equiv0$、$u_t=\Delta u$、而 $u(\cdot,0)\equiv0$ 的解，其增长快到 $e^{a|x|^2}$ 量级。

所以全空间热方程的正确提法是「在某个增长类里求唯一解」：

- **有界解**：唯一（用极值原理，即 §22.5 第 1 条的加权平均）；
- **$|u(x,t)|\le Ce^{a|x|^2}$、$a$ 足够小**：唯一（$a$ 的临界值与 $\frac{1}{4t}$ 同阶 —— 热核自身就是 $e^{-\frac{|x|^2}{4t}}$ 量级，增长再快就「跑得比扩散还快」）；
- **无任何增长限制**：不唯一，Tychonoff 反例。

这也是为什么 [[18-热方程·建模、边值条件与分离变量法]] 在有界区域里从没提过这回事：那里区域有界，$\sup|u|$ 自动有限。

---

## 22.7 非齐次 Cauchy 问题与 Duhamel 原理

(22.1) 是齐次方程。加上热源就是

$$
\begin{cases}\dfrac{\partial u}{\partial t}-\Delta u=f(x,t), & (x,t)\in\mathbb R^n\times(0,+\infty),\\[6pt] u(x,0)=\varphi(x), & x\in\mathbb R^n .\end{cases}\tag{22.8}
$$

结论：解为

$$
\boxed{\ u(x,t)=\underbrace{\frac{1}{(4\pi t)^{n/2}}\int_{\mathbb R^n}\varphi(y)e^{-\frac{|x-y|^2}{4t}}dy}_{=\,G_t*\varphi}+\underbrace{\int_0^t\!\!\int_{\mathbb R^n}\frac{f(y,\tau)}{(4\pi(t-\tau))^{n/2}}e^{-\frac{|x-y|^2}{4(t-\tau)}}dy\,d\tau}_{=\,\int_0^tG_{t-\tau}*f(\cdot,\tau)\,d\tau}\ }\tag{22.9}
$$

结构与在 [[20-Duhamel原理与齐次化]] 里完全一样 —— **齐次解 + 源项沿时间的叠加**，只是那里的「传播子」是 $e^{(t-\tau)A}$（特征展开），这里是热核 $G_{t-\tau}$。

**正确性（验证，而非推导）。** 只需验证第二项 $u_2:=\int_0^tG_{t-\tau}*f(\cdot,\tau)d\tau$ 满足齐次初值、零初值、右端 $f$。与 [[20-Duhamel原理与齐次化]] §20.3 同一套 Leibniz 法则：

$$
\partial_tu_2=\underbrace{\left.G_{t-\tau}*f(\cdot,\tau)\right|_{\tau=t}}_{=\,G_0*f(\cdot,t)}+\int_0^t\partial_tG_{t-\tau}*f(\cdot,\tau)\,d\tau .
$$

第一项是上限项：$\tau=t$ 时 $t-\tau=0$，而 $G_0$ 在 $t\to0^+$ 时是 $\delta_0$（§22.5 第 4 条），故 $G_0*f(\cdot,t)=f(\cdot,t)$ —— **这一项正是整个原理的唯一实质**，与 [[20-Duhamel原理与齐次化]] 订正第 1 条是同一个坑（那里错写成 $f(x,\tau)$ 把哑变量留在了结果里）。再由 $\partial_tG_t=\Delta G_t$（§22.5 第 3 条）与 $\Delta_x(G_{t-\tau}*f)=\Delta G_{t-\tau}*f$，

$$
\partial_tu_2=f(\cdot,t)+\int_0^t\Delta G_{t-\tau}*f(\cdot,\tau)\,d\tau=f(\cdot,t)+\Delta u_2 .
$$

初值：$u_2(x,0)=\int_0^0(\cdots)=0$ ✓。由线性性，$G_t*\varphi+u_2$ 满足 (22.8)。$\square$

> **要点（$G_0=\delta_0$ 这一句是全部关键）**　只要承认「$t\to0^+$ 时 $G_t$ 是近似单位、极限为 $\delta_0$」，Duhamel 原理就从有界区域原封不动搬到全空间。因为 $G_t$ 恰好是「从一点出发的热量分布」，$G_{t-\tau}*f(\cdot,\tau)$ 就是「$\tau$ 时刻注入的那一层热量，自由扩散 $t-\tau$ 时间后的样子」。**热核的半群性质 $G_t*G_s=G_{t+s}$（§22.5 第 2 条）保证了这种"接力式"解释自洽。**

> **补充（讲稿这里只写了一句话）**　讲稿的证明是「Duhamel 原理：解 $\partial_tw-\Delta w=0$、$w(x,0)=f(x,\tau)$，再积分」—— 只说了解法，没验。而且它没有交代两件事：① $w$ 里的 $x$ 与积分变量、与最终解的自变量之间的冲突（应像 [[20-Duhamel原理与齐次化]] 订正第 2 条那样改用 $w_\tau(x,s)$、$s=t-\tau$）；② 上限项为什么给 $f(x,t)$ 而不是 $f(x,\tau)$。本篇已按 $w_\tau(x,s)$ 的口径重写。

---

## 22.8 与前面内容的接口

**补充（有界区域的热核与听音辨鼓）**　把 (22.5) 与 19 篇的级数解对照。设 $\{(\lambda_k,u_k)\}$ 为 $-\Delta$（配 Dirichlet 条件）在 $\Omega$ 上的特征对，则

$$
p_\Omega(x,y;t):=\sum_{k=1}^{\infty}e^{-\lambda_kt}u_k(x)u_k(y)
$$

是 $\Omega$ 上的热核（$u(x,t)=\int_\Omega p_\Omega(x,y;t)\varphi(y)dy$ 解初边值问题）。它与 $G_t$ 的关系：$p_\Omega$ 是 $G_t$ 在 $\Omega\times\Omega$ 上「用 $\Omega$ 的特征函数做的谱投影」，$t\to0^+$ 时 $p_\Omega(x,y;t)\to\delta(x-y)$，与 $G_t$ 的差别只在边界的镜像修正。

于是 Laplace 谱 $\{\lambda_k\}$ 完全由 $p_\Omega$ 决定，而 $p_\Omega$ 由「$\Omega$ 上热传导的观测」决定。1966 年 Kac 提出「**能否听出鼓的形状**」：$\mathrm{tr}\,p_\Omega(t)=\sum_ke^{-\lambda_kt}$ 已知时能否定出 $\Omega$？在 $\mathbb R^2$ 上 $t\to0^+$ 的渐近展开给出

$$
\sum_ke^{-\lambda_kt}=\frac{|\Omega|}{4\pi t}\mp\frac{|\partial\Omega|}{8\sqrt{\pi t}}+O(1),
$$

即只能定出**面积与周长**；1992 年 Gordon–Webb–Wolpert 造出等谱而不等距的平面区域，答案是否定的。（细节见 [[14-谱理论·鼓面的数学建模与Laplace特征值]] §14.6。）

**补充（与基本解、Green 函数的关系）**　(22.7) 对 $t$ 从 $0$ 到 $\infty$ 积分，得

$$
\int_0^{\infty}G_t(x)\,dt=\frac{\Gamma\left(\frac n2-1\right)}{4\pi^{n/2}}|x|^{2-n}\qquad(n\ge3),
$$

右端正是 $-\Delta$ 在 $\mathbb R^n$ 的基本解 $\Gamma_n(x)$（[[07-基本解与广义函数]]）。即

$$
(-\Delta)^{-1}=\int_0^{\infty}e^{t\Delta}\,dt ,
$$

「椭圆算子的逆 = 抛物半群沿时间的积分」。同理，有界区域上 $\int_0^\infty p_\Omega(x,y;t)dt=G_\Omega(x,y)$ 就是 Dirichlet Green 函数（[[09-Green函数的性质]]）—— 两条路线（Green 表示公式、热核积分）给出同一个东西。

---

## 订正

1. **热核计算的第一行漏了 $\dfrac{1}{(2\pi)^n}$（第八周版本）。** 讲稿写 $\mathcal F^{-1}\left[e^{-|\xi|^2t}\right](x)=\int_{\mathbb R^n}e^{-|\xi|^2t+\sqrt{-1}\langle\xi,x\rangle}d\xi$，按本课约定（[[21-Fourier变换与速降函数空间]] 定义 21.3.1）逆变换必须带 $\frac{1}{(2\pi)^n}$。这不是记号的洁癖：漏掉它，再配上正确的 Gauss 积分 $(\pi/t)^{n/2}$，得到的是 $\frac{1}{2^n\pi^n}\frac{\pi^{n/2}}{t^{n/2}}$，与真值 $\frac{1}{2^n\pi^{n/2}t^{n/2}}$ 差一个因子 $\pi^{-n/2}$，最终 $\int G_t\ne1$。本篇 §22.4 已补。**第九周重述这一段时已把这个因子补上**，两版对照即可看出是排印疏漏而非概念错误。

2. **第二行把 $\dfrac{|x|^2}{4t}$ 算了两次 —— 这是本篇最要紧的一处。** 讲稿写
$\left[\int e^{-t\left(|\xi|^2-\frac{\sqrt{-1}\langle\xi,x\rangle}{t}+\frac{|x|^2}{4t^2}\right)}d\xi\right]e^{-\frac{|x|^2}{4t}}$，
括号里那个 $+\frac{|x|^2}{4t^2}$ 乘上 $-t$ 就是 $-\frac{|x|^2}{4t}$，已经产生了外面的那个因子；再乘一次 $e^{-\frac{|x|^2}{4t}}$，总指数变成 $-\frac{|x|^2}{2t}$。**判据**：把它算完得 $G_t=\frac{1}{(2\pi)^n}(\pi/t)^{n/2}e^{-\frac{|x|^2}{2t}}$，于是 $\int_{\mathbb R^n}G_t=2^{-n}$，$n=1$ 时是 $\frac12$ —— 热量凭空少了一半。正确的中间步只有两种写法：要么把配方的结果整块写成
$-t|\xi|^2+\sqrt{-1}\langle\xi,x\rangle=-t\left|\xi-\frac{\sqrt{-1}x}{2t}\right|^2-\frac{|x|^2}{4t}$
再一次性提出 $e^{-\frac{|x|^2}{4t}}$（本篇 §22.4.1 的 (22.6)）；要么保留括号内的 $+\frac{|x|^2}{4t^2}$、但外面只写 $e^{-\frac{|x|^2}{4t}}$ 而**不**先把它提出来。讲稿把两种写法混在了一行里。**第九周重述这一段时已改对**：新版括号内是 $\sum_j\xi_j^2-\frac{\sqrt{-1}}{t}\sum_jx_j\xi_j+\left(\frac{\sqrt{-1}}{2t}\right)^2\sum_jx_j^2$、外面配 $e^{t\left(\frac{\sqrt{-1}}{2t}\right)^2\sum_jx_j^2}=e^{-\frac{|x|^2}{4t}}$，与本节 (22.6) 等价。**所以第 1、2 两条只针对第八周那份讲稿**，两版对照阅读时不必以为讲稿前后自相矛盾 —— 它是把同一段推导重写了一遍。

3. **「也假设 $u(\cdot,t)\in\mathcal S(\mathbb R^n)$」是循环论证。** 讲稿先假设解在 $\mathcal S$ 里（才能对 $u$ 作 $\mathcal F$），再由方程解出 $\widehat u$，最后声称这就是解 —— 但「解存在且属于 $\mathcal S$」本身从未验证。正确顺序是**先构造后验证**：令 $u:=G_t*\varphi$（$\varphi\in\mathcal S$ 保证 $G_t*\varphi\in\mathcal S$，因为 $G_t\in\mathcal S$ 且 $\mathcal S$ 对卷积封闭），再逐条验证它满足 (22.1)（用 §22.5 第 3 条与命题 22.6.2）。这样 $\mathcal S$ 只是「让所有积分与分部积分合法」的技术假设，而不是关于解的先验断言。本篇 §22.4–22.6 按这个顺序写。

4. **$\mathcal F\left[\dfrac{\partial u}{\partial t}\right]=\dfrac{\partial}{\partial t}\mathcal F[u]$ 未给理由。** 需要 $\partial_tu$ 关于 $x$ 绝对可积并且积分号下求导合法。在 $\mathcal S$ 框架下这可以由 $\mathcal F$ 的显式表达式直接看出：$\widehat u(\xi,t)=\int u(x,t)e^{-\sqrt{-1}\langle x,\xi\rangle}dx$，若 $u,\partial_tu$ 都速降（对 $x$ 一致），则 $|\partial_tu|\le C_k(1+|x|)^{-k}$，控制收敛即给出可交换。讲稿只写「则」字带过。

5. **卷积那一步的记号会让人误解。** 讲稿写 $u(x,t)=\mathcal F^{-1}\left[e^{-|\xi|^2t}\right](x)*\varphi(x)$，把「$\mathcal F^{-1}[\cdots]$ 在 $x$ 处的值」与 $\varphi(x)$ 作卷积 —— 字面上是两个数相乘，而且卷积的两个因子都写了自变量 $x$。应为 $(G_t*\varphi)(x)$，先把 $G_t:=\mathcal F^{-1}\left[e^{-|\xi|^2t}\right]$ 命名，再写卷积。这个命名还关系到标题里的「热核」二字：讲稿从头到尾没给热核定名，标题与内容脱节。本篇 §22.3 补了定义。

6. **「$\varphi\in\mathcal S$」比实际需要强得多。** 由命题 22.6.1 与 22.6.2，$\varphi\in L^\infty$ 就够（$u$ 在 $t>0$ 时 $C^\infty$、$t\to0^+$ 时在连续点收敛）。$\mathcal S$ 的价值在于让 $\mathcal F$ 与卷积在同一个空间里闭合，属于技术便利，不是本质限制。讲稿把它写成问题的前提，会让读者以为光滑性是解存在的必要条件 —— 恰恰相反，热方程最能容忍坏的初值。

---

系列：[[21-Fourier变换与速降函数空间]] · [[18-热方程·建模、边值条件与分离变量法]] · [[19-分离变量解的正则性·三个收敛命题]]
