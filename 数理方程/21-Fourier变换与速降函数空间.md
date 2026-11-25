# 21 Fourier 变换与速降函数空间

> 系列：[[20-Duhamel原理与齐次化]] · [[22-Cauchy问题与热核]] · [[15-Sturm-Liouville定理与一维特征值问题]]

到此为止，一切「展开」都建立在 $(0,L)$ 上：特征值离散、特征函数是三角函数（[[15-Sturm-Liouville定理与一维特征值问题]]）。这一篇把 $L\to\infty$ 取掉，看会发生什么 —— 结果是**离散谱变成连续谱**，级数变成积分，而「特征函数」变成 $e^{\sqrt{-1}\langle\xi,x\rangle}$。

讲稿的过渡只有一句话：「随着 $L$ 增大，$\lambda_k$ 与 $\lambda_{k+1}$ 之间的间距不断变小，最终趋向于 $0$」。这句话值得展开，因为它是整门课里唯一一次「从可数到不可数」的跳变。

---

## 21.1 从离散谱到连续谱

在 $(0,L)$ 上（Dirichlet），

$$
\lambda_k=\left(\frac{k\pi}{L}\right)^2,\qquad k=1,2,\dots;\qquad \lambda_{k+1}-\lambda_k=\frac{\pi^2}{L^2}\left(2k+1\right).
$$

固定 $k$（看谱的「低端」），相邻间距 $\sim\dfrac{2k\pi^2}{L^2}\to0$ 当 $L\to\infty$。**谱越来越密，极限情形下成为连续统。**

**极限情形的「特征函数」**是什么？在 $\mathbb R^n$ 上试着找 $-\Delta$ 的特征函数：对于参数 $\xi\in\mathbb R^n$ 令

$$
u_\xi(x):=e^{\sqrt{-1}\langle\xi,x\rangle},
$$

直接算（每求一次二阶导就掉下一个 $(\sqrt{-1}\xi_j)^2=-\xi_j^2$）：

$$
-\Delta e^{\sqrt{-1}\langle\xi,x\rangle}=\left(\sum_{j=1}^{n}\xi_j^2\right)e^{\sqrt{-1}\langle\xi,x\rangle}=|\xi|^2e^{\sqrt{-1}\langle\xi,x\rangle}.
$$

所以 $u_\xi$ 是特征值 $|\xi|^2$ 的特征函数，**谱是 $[0,+\infty)$ 整个连续统**（每个 $|\xi|^2$ 对应无穷多个 $\xi$，即特征空间无穷维 —— 离散情形的「重数有限」在这里彻底失效）。

**困难出在边值条件上。** 有界区域的情形靠「$\partial\Omega$ 上取零」把谱筛成离散；$\mathbb R^n$ 没有边界，无处施加条件。而 $u_\xi$ 本身又不可积（$|u_\xi|\equiv1$）。讲稿的处理是：**不给 $u$ 加条件，改给 $f$ 加条件** —— 要求待变换的函数属于速降函数空间 $\mathcal S(\mathbb R^n)$，它衰减得足够快，保证一切积分收敛、一切分部积分合法。

> **要点**　这一步是「用衰减换边值」：有界区域用**边界条件**制造离散谱；全空间用**衰减条件**让积分有意义。两者不是同一件事，但作用类似 —— 都是为了让「特征展开」收敛。

---

## 21.2 速降函数空间 $\mathcal S(\mathbb R^n)$

**定义 21.2.1**　称 $\mathcal S(\mathbb R^n)$ 为**速降函数空间**（或 Schwartz 空间）：

$$
\mathcal S(\mathbb R^n):=\left\{\varphi\in C^\infty(\mathbb R^n)\ \Big|\ \lim_{|x|\to+\infty}(1+|x|^k)|D^\alpha\varphi(x)|=0,\quad \forall k\in\mathbb N_0,\ \forall\alpha\in\mathbb N_0^n\right\},
$$

其中 $\alpha=(\alpha_1,\dots,\alpha_n)$、$D^\alpha\varphi=\partial_{x_1}^{\alpha_1}\cdots\partial_{x_n}^{\alpha_n}\varphi$、$|\alpha|=\sum_j\alpha_j$。

**例**

- $e^{-|x|^2}\in\mathcal S(\mathbb R^n)$：任意次求导后只多出多项式因子，$e^{-|x|^2}$ 把它们全压住；
- $C_c^\infty(\mathbb R^n)\subseteq\mathcal S(\mathbb R^n)$：紧支撑函数在无穷远处直接为 $0$。

> **要点（这个空间为什么是「对的」）**　$\mathcal S$ 的三条性质恰好够用：
> 1. **对求导封闭**：$D^\alpha\varphi\in\mathcal S$；
> 2. **对乘多项式封闭**：$x^\beta\varphi\in\mathcal S$；
> 3. **$\mathcal S\subseteq L^1(\mathbb R^n)$**：取 $k=n+1$，$(1+|x|)^{n+1}|\varphi|$ 有界，故 $\int|\varphi|\lesssim\int(1+|x|)^{-(n+1)}<\infty$。
>
> 第 3 条让 Fourier 积分绝对收敛，第 1、2 条让「分部积分把导数搬到另一边」永远不用管边界项（无穷远处全为零）。**这就是为什么 Fourier 变换理论要在 $\mathcal S$ 上做**：在 $L^2$ 上做，分部积分立刻要处理边界。

---

## 21.3 Fourier 变换与逆变换

**定义 21.3.1**　对 $f\in\mathcal S(\mathbb R^n)$，定义 **Fourier 变换**

$$
\mathcal F[f](\xi):=\int_{\mathbb R^n}f(x)e^{-\sqrt{-1}\langle x,\xi\rangle}\,dx ,
$$

及 **Fourier 逆变换**

$$
\mathcal F^{-1}[g](x):=\frac{1}{(2\pi)^n}\int_{\mathbb R^n}g(\xi)e^{\sqrt{-1}\langle\xi,x\rangle}\,d\xi .
$$

两者互逆（$\mathcal F^{-1}\mathcal F=\mathrm{id}$，$\mathcal F\mathcal F^{-1}=\mathrm{id}$），这就是 §21.1 里那个「连续版特征展开」：

$$
f(x)=\frac{1}{(2\pi)^n}\int_{\mathbb R^n}\widehat f(\xi)\,u_\xi(x)\,d\xi,\qquad \widehat f:=\mathcal F[f].
$$

把它与有界区域的展开 $f=\sum_k(f,u_k)u_k$ 对比：**求和变积分、$\delta_{k\ell}$ 变 $\delta(\xi-\eta)$、系数 $(f,u_k)$ 变成 $\widehat f(\xi)$**。

> **补充（约定的自由度，以及为什么本课选这一种）**　常见的三种约定：
>
> | 约定 | 正变换 | 逆变换 | Plancherel |
> |---|---|---|---|
> | 本课 | $\int fe^{-\sqrt{-1}\langle x,\xi\rangle}dx$ | $\frac{1}{(2\pi)^n}\int ge^{\sqrt{-1}\langle\xi,x\rangle}d\xi$ | $\int f\bar g=\frac{1}{(2\pi)^n}\int\widehat f\overline{\widehat g}$ |
> | 对称 | $(2\pi)^{-n/2}\int fe^{-\sqrt{-1}\langle x,\xi\rangle}dx$ | $(2\pi)^{-n/2}\int ge^{\sqrt{-1}\langle\xi,x\rangle}d\xi$ | $\int f\bar g=\int\widehat f\overline{\widehat g}$ |
> | 无 $2\pi$ | $\int fe^{-2\pi\sqrt{-1}\langle x,\xi\rangle}dx$ | $\int ge^{2\pi\sqrt{-1}\langle\xi,x\rangle}d\xi$ | $\int f\bar g=\int\widehat f\overline{\widehat g}$ |
>
> 本课选第一种（$2\pi$ 全给逆变换），好处是与 $\mathcal F[\partial^\alpha f]=(\sqrt{-1}\xi)^\alpha\widehat f$ 配合时最干净；代价是 Plancherel 里多一个 $(2\pi)^{-n}$。**换约定时只需要盯住两件事：$2\pi$ 的分布、以及指数上的 $\pm$ 号。** 后面 [[22-Cauchy问题与热核]] 的热核常数就取决于这一选择。

### 21.3.1 互逆性的证明：Gauss 积分与扰动法

上面说「两者互逆」，但到此为止只是断言。讲稿第九周给了完整证明，思路是**先算带参数 $\delta$ 的 Gauss 积分，再用 $\delta\to0$ 的扰动把不绝对收敛的积分合法化**。这套路子在 [[22-Cauchy问题与热核]] 里还会以「$G_t$ 是近似单位」的形式再出现一次，值得记牢。

**引理 21.3.2（Gauss 积分公式）**　对 $\varepsilon>0$、$x\in\mathbb R^n$，

$$
G_\varepsilon(x):=\int_{\mathbb R^n}e^{-\frac{\varepsilon}{2}|\xi|^2+\sqrt{-1}\langle x,\xi\rangle}d\xi=\left(\frac{2\pi}{\varepsilon}\right)^{n/2}e^{-\frac{|x|^2}{2\varepsilon}}.\tag{21.5}
$$

**证明**　被积函数关于 $\xi$ 是「指数衰减 $\times$ 模为 $1$ 的振荡」，故 $G_\varepsilon$ 可微且导数可搬进积分号。对 $x_j$ 求导，然后把 $\sqrt{-1}\xi_j e^{-\frac\varepsilon2|\xi|^2}$ 写成对 $\xi_j$ 的导数、分部积分（边界项因指数衰减消失）：

$$
\partial_{x_j}G_\varepsilon=\int_{\mathbb R^n}\sqrt{-1}\xi_j\,e^{-\frac\varepsilon2|\xi|^2+\sqrt{-1}\langle x,\xi\rangle}d\xi=-\frac{\sqrt{-1}}{\varepsilon}\int_{\mathbb R^n}\frac{\partial}{\partial\xi_j}\left(e^{-\frac\varepsilon2|\xi|^2}\right)e^{\sqrt{-1}\langle x,\xi\rangle}d\xi=-\frac{x_j}{\varepsilon}G_\varepsilon .
$$

（末步：分部积分把 $\partial_{\xi_j}$ 搬到 $e^{\sqrt{-1}\langle x,\xi\rangle}$ 上，得 $-\sqrt{-1}x_j$，而 $(-\sqrt{-1})(-\sqrt{-1})=1$。）于是 $G_\varepsilon$ 满足 $\nabla G_\varepsilon=-\frac{x}{\varepsilon}G_\varepsilon$，解得 $G_\varepsilon(x)=G_\varepsilon(0)e^{-\frac{|x|^2}{2\varepsilon}}$。再算常数：

$$
G_\varepsilon(0)=\int_{\mathbb R^n}e^{-\frac\varepsilon2|\xi|^2}d\xi=\left(\sqrt{\frac{2\pi}{\varepsilon}}\right)^n=\left(\frac{2\pi}{\varepsilon}\right)^{n/2}. \qquad\square
$$

**第一步：先把 $\widehat f$ 弄成 $L^1$ 的。** 反演公式的右端是 $\frac{1}{(2\pi)^n}\int\widehat f(\xi)e^{\sqrt{-1}\langle\xi,y\rangle}d\xi$，要让它有意义，$\widehat f$ 必须绝对可积。取整数 $N>n/2$，由 $\mathcal F[(1-\Delta)^Nf](\xi)=(1+|\xi|^2)^N\widehat f(\xi)$（命题 21.4.1 的反复应用），

$$
\left|\widehat f(\xi)\right|\le\left\lVert(1-\Delta)^Nf\right\rVert_{L^1(\mathbb R^n)}\left(1+|\xi|^2\right)^{-N},
$$

而 $\int(1+|\xi|^2)^{-N}d\xi<+\infty$（$2N>n$），故 $\widehat f\in L^1(\mathbb R^n)$。

**第二步：扰动。** 直接算 $\mathcal F^{-1}\mathcal F$ 会碰到 $\iint|f(x)|\,dx\,d\xi=+\infty$ 这个死结（$e^{\sqrt{-1}\langle y-x,\xi\rangle}$ 模为 $1$，对 $\xi$ 不可积），**Fubini 不成立**。所以插入一个衰减因子 $e^{-\frac\varepsilon2|\xi|^2}$：

$$
I_\varepsilon(y):=\frac{1}{(2\pi)^n}\int_{\mathbb R^n}\left[\int_{\mathbb R^n}e^{\sqrt{-1}\langle y-x,\xi\rangle-\frac{\varepsilon}{2}|\xi|^2}d\xi\right]f(x)\,dx ,\qquad \varepsilon>0 .\tag{21.6}
$$

现在内层积分就是引理 21.3.2（令 $\xi\mapsto-\xi$，注意 $\langle y-x,-\xi\rangle=\langle x-y,\xi\rangle$，结果不变）：

$$
I_\varepsilon(y)=\frac{1}{(2\pi)^n}\int_{\mathbb R^n}\left(\frac{2\pi}{\varepsilon}\right)^{n/2}e^{-\frac{|y-x|^2}{2\varepsilon}}f(x)\,dx=\frac{1}{(2\pi\varepsilon)^{n/2}}\int_{\mathbb R^n}e^{-\frac{|x|^2}{2\varepsilon}}f(y-x)\,dx ,
$$

末步用了 $\left|y-x\right|$ 的对称性与平移（把 $x$ 换成 $y-x$ 不改变 $\mathbb R^n$ 上的积分）。再作换元 $x=\sqrt{2\varepsilon}\,\zeta$：

$$
I_\varepsilon(y)=\frac{1}{\pi^{n/2}}\int_{\mathbb R^n}e^{-|\zeta|^2}f\left(y-\sqrt{2\varepsilon}\,\zeta\right)d\zeta .\tag{21.7}
$$

**第三步：$\varepsilon\to0^+$ 的两种算法。**

*算法一（看 (21.7)）。* $\lvert f\rvert\le M$ 且 $f$ 连续，被积函数被可积函数 $Me^{-|\zeta|^2}$ 控制，由控制收敛定理

$$
\lim_{\varepsilon\to0^+}I_\varepsilon(y)=\frac{f(y)}{\pi^{n/2}}\int_{\mathbb R^n}e^{-|\zeta|^2}d\zeta=f(y) ,
$$

最后一步用了 $\int_{\mathbb R^n}e^{-|\zeta|^2}d\zeta=\pi^{n/2}$（即 (21.5) 取 $\varepsilon=2$、$x=0$）。

*算法二（回到 (21.6)，先对 $x$ 积分）。* 内层对 $x$ 的积分正是 $\widehat f(\xi)$，故

$$
I_\varepsilon(y)=\frac{1}{(2\pi)^n}\int_{\mathbb R^n}\widehat f(\xi)e^{\sqrt{-1}\langle\xi,y\rangle-\frac{\varepsilon}{2}|\xi|^2}d\xi .
$$

由第一步 $\widehat f\in L^1$，且被积函数被 $|\widehat f|$ 控制、逐点趋于 $\widehat f(\xi)e^{\sqrt{-1}\langle\xi,y\rangle}$，控制收敛定理给

$$
\lim_{\varepsilon\to0^+}I_\varepsilon(y)=\frac{1}{(2\pi)^n}\int_{\mathbb R^n}\widehat f(\xi)e^{\sqrt{-1}\langle\xi,y\rangle}d\xi=\mathcal F^{-1}\left[\widehat f\right](y) .
$$

**结论。** 两种算法给出同一极限，故 $\mathcal F^{-1}\mathcal F[f](y)=f(y)$。$\square$

> **要点（扰动法为什么必须分两步看同一个量）**　这是全篇最"绕"的一点：$I_\varepsilon$ **一个**对象，被算成**两个**表达式，各自取极限，极限相等 —— 于是把"未知的 $\mathcal F^{-1}\mathcal F$"和"已知的 $f$"接上了。**衰减因子 $e^{-\frac\varepsilon2|\xi|^2}$ 的作用是让两个方向都可积**：对 $\xi$ 方向，它把振荡积分变成 Gauss 积分（引理 21.3.2 可用）；对 $x$ 方向，它把 $f$ 的积分变成近似单位卷积（$f$ 的连续性可用）。缺一个方向就断链。
>
> 顺带：$\frac{1}{\pi^{n/2}}e^{-|\zeta|^2}$ 就是方差 $1/2$ 的 Gauss 密度，所以 (21.7) 说的是「$I_\varepsilon$ 是 $f$ 的一个近似单位卷积」—— 这正是 [[22-Cauchy问题与热核]] §22.5 第 4 条的同一件事，只是那里用的是 $H_t$、参数是 $t$ 而不是 $\varepsilon$。**Fourier 反演与「热核还原初值」是同一个极限。**

---

## 21.4 基本性质

### 21.4.1 线性

任取 $f_1,f_2\in\mathcal S$、$a,b\in\mathbb C$，由积分的线性性，

$$
\mathcal F[af_1+bf_2](\xi)=a\mathcal F[f_1](\xi)+b\mathcal F[f_2](\xi).
$$

即 $\mathcal F:\mathcal S(\mathbb R^n)\to\mathcal S(\mathbb R^n)$ 是线性映射（为什么会落在 $\mathcal S$ 里见 §21.4.4）。

### 21.4.2 求导变乘法

**命题 21.4.1**　对 $f\in\mathcal S(\mathbb R^n)$，

$$
\mathcal F[\partial^\alpha f](\xi)=(\sqrt{-1}\xi)^\alpha\mathcal F[f](\xi).
$$

**证明**　先做一阶、第一个变量：

$$
\int_{\mathbb R^n}\partial_{x_1}f(x)e^{-\sqrt{-1}\langle x,\xi\rangle}dx=\int_{\mathbb R^n}\partial_{x_1}\left(f(x)e^{-\sqrt{-1}\langle x,\xi\rangle}\right)dx+\sqrt{-1}\xi_1\int_{\mathbb R^n}f(x)e^{-\sqrt{-1}\langle x,\xi\rangle}dx ,
$$

（因为 $\partial_{x_1}\left(fe^{-\sqrt{-1}\langle x,\xi\rangle}\right)=\partial_{x_1}f\,e^{-\sqrt{-1}\langle x,\xi\rangle}-\sqrt{-1}\xi_1fe^{-\sqrt{-1}\langle x,\xi\rangle}$，移项即得。）第一项是 $x_1$ 方向的全导数积分：

$$
\int_{\mathbb R^n}\partial_{x_1}(\cdots)dx=\int_{\mathbb R^{n-1}}\left[\int_{-\infty}^{+\infty}\frac{\partial}{\partial x_1}\left(fe^{-\sqrt{-1}\langle x,\xi\rangle}\right)dx_1\right]dx'=0,
$$

因为 $f\in\mathcal S$ 使 $fe^{-\sqrt{-1}\langle x,\xi\rangle}$ 在 $x_1\to\pm\infty$ 时趋于 $0$（其模长就是 $|f(x)|$，而 $f$ 速降）。故 $\mathcal F[\partial_{x_1}f]=\sqrt{-1}\xi_1\mathcal F[f]$。一般 $\alpha$ 由归纳（$\partial^\alpha f\in\mathcal S$，可以逐次用上面的一阶公式）即得。$\square$

### 21.4.3 乘 $x$ 变求导

**命题 21.4.2**　对 $f\in\mathcal S(\mathbb R^n)$，

$$
\mathcal F[x^\alpha f(x)](\xi)=(\sqrt{-1})^\alpha\partial_\xi^\alpha\mathcal F[f](\xi).
$$

**证明**　对 $\mathcal F[f](\xi)=\int f(x)e^{-\sqrt{-1}\langle x,\xi\rangle}dx$ 关于 $\xi_j$ 求导，可把导数搬进积分号（$f$ 速降保证一致收敛）：$\partial_{\xi_j}e^{-\sqrt{-1}\langle x,\xi\rangle}=-\sqrt{-1}x_je^{-\sqrt{-1}\langle x,\xi\rangle}$，故

$$
\partial_{\xi_j}\mathcal F[f](\xi)=-\sqrt{-1}\int x_jf(x)e^{-\sqrt{-1}\langle x,\xi\rangle}dx=-\sqrt{-1}\,\mathcal F[x_jf](\xi),
$$

即 $\mathcal F[x_jf]=\sqrt{-1}\partial_{\xi_j}\mathcal F[f]$。反复用即得。$\square$

> **要点**　21.4.2 与 21.4.3 是一对：**求导 ↔ 乘坐标**在变换下互换。这是 Fourier 变换能被用来解常系数 PDE 的全部原因 —— 它把 $\Delta$ 变成乘以 $|\xi|^2$，把微分方程变成代数方程。后面 [[22-Cauchy问题与热核]] 就是直接受益者。

### 21.4.4 $\mathcal F$ 把 $\mathcal S$ 映到 $\mathcal S$

**定理 21.4.1**　若 $f\in\mathcal S(\mathbb R^n)$，则 $\widehat f:=\mathcal F[f]\in\mathcal S(\mathbb R^n)$。

**证明**　要证对任意 $\alpha,\beta\in\mathbb N_0^n$ 都有 $\sup_{\xi}|\xi^\alpha\partial_\xi^\beta\widehat f(\xi)|<+\infty$。由 21.4.2 与 21.4.3：

$$
(\sqrt{-1}\xi)^\alpha\partial_\xi^\beta\widehat f(\xi)=(\sqrt{-1}\xi)^\alpha\mathcal F\left[(-\sqrt{-1}x)^\beta f\right](\xi)=(-\sqrt{-1})^{|\beta|}\mathcal F\left[\partial_x^\alpha(x^\beta f)\right](\xi).
$$

（第一步用 21.4.2，第二步用 21.4.2 的反向：$(\sqrt{-1}\xi)^\alpha\mathcal F[g]=\mathcal F[\partial^\alpha g]$。）而 $\partial_x^\alpha(x^\beta f)\in\mathcal S\subseteq L^1(\mathbb R^n)$（§21.2 要点第 3 条），故

$$
\left|\xi^\alpha\partial_\xi^\beta\widehat f(\xi)\right|=\left|\mathcal F\left[\partial_x^\alpha(x^\beta f)\right](\xi)\right|\le\left\lVert\partial_x^\alpha(x^\beta f)\right\rVert_{L^1(\mathbb R^n)}<+\infty ,
$$

与 $\xi$ 无关。$\square$

> **补充（定理 21.4.1 只证了一半）**　讲稿的标题说「$\mathcal F$ 与 $\mathcal F^{-1}$ 给出 $\mathcal S$ 到 $\mathcal S$ 的**同构**」，但证明只做了「$\mathcal F$ 把 $\mathcal S$ 映进 $\mathcal S$」。另一半「$\mathcal F^{-1}$ 也把 $\mathcal S$ 映进 $\mathcal S$」可以由同一结论推出：把 $\xi$ 换成 $-\xi$ 即得

$$
\mathcal F^{-1}[g](x)=\frac{1}{(2\pi)^n}\mathcal F[g](-x),
$$

而 $\mathcal S$ 对「反射 $x\mapsto-x$」封闭，故 $\mathcal F^{-1}$ 的表现与 $\mathcal F$ 只差一个常数与反射。再由 $\mathcal F^{-1}\mathcal F=\mathrm{id}$ 得双射。

### 21.4.5 卷积定理

**定义 21.4.3**　对 $f_1,f_2\in\mathcal S$，称

$$
(f_1*f_2)(x):=\int_{\mathbb R^n}f_1(x-y)f_2(y)\,dy=\int_{\mathbb R^n}f_1(y)f_2(x-y)\,dy
$$

为 $f_1$ 与 $f_2$ 的**卷积**（两个表达式相等由换元 $y\mapsto x-y$ 得到）。

**定理 21.4.2**　若 $f_1,f_2\in\mathcal S(\mathbb R^n)$，则

$$
\mathcal F[f_1*f_2](\xi)=\mathcal F[f_1](\xi)\cdot\mathcal F[f_2](\xi),\qquad
\mathcal F^{-1}[f_1\cdot f_2](x)=\left(\mathcal F^{-1}[f_1]*\mathcal F^{-1}[f_2]\right)(x).
$$

**证明**　直接算，把指数拆开 $e^{-\sqrt{-1}\langle x,\xi\rangle}=e^{-\sqrt{-1}\langle x-y,\xi\rangle}e^{-\sqrt{-1}\langle y,\xi\rangle}$：

$$
\begin{aligned}
\mathcal F[f_1*f_2](\xi)&=\int_{\mathbb R^n}\left(\int_{\mathbb R^n}f_1(x-y)f_2(y)\,dy\right)e^{-\sqrt{-1}\langle x,\xi\rangle}dx\\
&=\int_{\mathbb R^n}\left(\int_{\mathbb R^n}f_1(x-y)e^{-\sqrt{-1}\langle x-y,\xi\rangle}dx\right)f_2(y)e^{-\sqrt{-1}\langle y,\xi\rangle}dy\\
&=\left(\int_{\mathbb R^n}f_1(z)e^{-\sqrt{-1}\langle z,\xi\rangle}dz\right)\left(\int_{\mathbb R^n}f_2(y)e^{-\sqrt{-1}\langle y,\xi\rangle}dy\right).
\end{aligned}
$$

第二个等号是 Fubini（$f_1,f_2$ 速降，绝对可积，合法），第三个等号是换元 $z=x-y$（内层积分与 $y$ 无关）。第二条由 $\mathcal F^{-1}$ 的线性与第一条等价。$\square$

> **要点**　卷积定理是「**卷积 ↔ 乘积**」的字典：微分方程里的非齐次卷积表示，在 Fourier 侧变成普通的乘积，可以逐点相除。[[22-Cauchy问题与热核]] 里解 $u=\mathcal F^{-1}\left[e^{-|\xi|^2t}\right]*\varphi$ 就是这条定理的直接应用 —— 而「求 $\mathcal F^{-1}\left[e^{-|\xi|^2t}\right]$」正是热核。

---

## 21.5 Plancherel 等式

**定理 21.5.1（Plancherel）**　若 $f,g\in\mathcal S(\mathbb R^n)$，则

$$
\int_{\mathbb R^n}f\bar g\,dx=\frac{1}{(2\pi)^n}\int_{\mathbb R^n}\mathcal F[f](\xi)\overline{\mathcal F[g](\xi)}\,d\xi .
$$

特别地 $\lVert\mathcal F[f]\rVert_{L^2}=(2\pi)^{n/2}\lVert f\rVert_{L^2}$。

**引理 21.5.2**　对 $f,h\in\mathcal S(\mathbb R^n)$，

$$
\int_{\mathbb R^n}\mathcal F[f](\xi)\,\overline{h(\xi)}\,d\xi=(2\pi)^n\int_{\mathbb R^n}f(x)\,\overline{\mathcal F^{-1}[h]}(x)\,dx .
$$

**证明（引理）**　把 $\mathcal F[f]$ 代入并用 Fubini，再取共轭：

$$
\int_{\mathbb R^n}\mathcal F[f](\xi)\overline{h(\xi)}\,d\xi=\int_{\mathbb R^n}f(x)\left[\int_{\mathbb R^n}\overline{h(\xi)}\,e^{-\sqrt{-1}\langle x,\xi\rangle}d\xi\right]dx .
$$

而

$$
\overline{\mathcal F^{-1}[h]}(x)=\overline{\frac{1}{(2\pi)^n}\int_{\mathbb R^n}h(\xi)e^{\sqrt{-1}\langle\xi,x\rangle}d\xi}=\frac{1}{(2\pi)^n}\int_{\mathbb R^n}\overline{h(\xi)}\,e^{-\sqrt{-1}\langle x,\xi\rangle}d\xi ,
$$

（用 $\langle\xi,x\rangle=\langle x,\xi\rangle$），故方括号 $=(2\pi)^n\overline{\mathcal F^{-1}[h]}(x)$。$\square$

**证明（定理）**　在引理中取 $h=\mathcal F[g]$，由 $\mathcal F^{-1}\mathcal F=\mathrm{id}$ 得 $\overline{\mathcal F^{-1}[\mathcal F[g]]}=\bar g$，于是

$$
\int_{\mathbb R^n}\mathcal F[f]\,\overline{\mathcal F[g]}\,d\xi=(2\pi)^n\int_{\mathbb R^n}f\bar g\,dx .
$$

两边除以 $(2\pi)^n$ 即得。取 $g=f$ 得范数等式。$\square$

> **要点（Plancherel 说什么）**　它说 $\mathcal F$ 在 $L^2$ 上是「几乎等距」的：差一个常数 $(2\pi)^{n/2}$。因为 $\mathcal S$ 在 $L^2$ 中稠密，还可以把 $\mathcal F$ **唯一地延拓**成 $L^2(\mathbb R^n)$ 上的有界算子（这叫 Plancherel 定理的延拓形式）—— 于是平方可积函数的 Fourier 变换总有意义，哪怕积分 $\int f e^{-\sqrt{-1}\langle x,\xi\rangle}dx$ 本身不绝对收敛。**这一步是把 Fourier 变换从「$\mathcal S$ 上的工具」升级成「$L^2$ 上的算子」的关键。**
>
> 若换成对称约定（§21.3 的表第二行），常数 $(2\pi)^{n/2}$ 消失，$L^2$ 上就是严格等距 —— 这是量子力学文献偏爱对称约定的原因。

---

## 21.6 与前面内容的接口

**补充（$\mathcal S$ 与广义函数）**　$\mathcal S$ 的真正价值在于它的对偶空间 $\mathcal S'$（**缓增广义函数**）：一切「不太坏」的广义函数（$\delta$、多项式、$e^{\sqrt{-1}\langle\xi,x\rangle}$）都在 $\mathcal S'$ 里，而 $\mathcal F$ 可以按「转置」定义到 $\mathcal S'$ 上：$\langle\mathcal F[T],\varphi\rangle:=\langle T,\mathcal F[\varphi]\rangle$。于是 $\mathcal F[\delta_0]=1$、$\mathcal F[1]=(2\pi)^n\delta_0$ 都有了严格含义。这条线接到 [[07-基本解与广义函数]]：那里用「$-\Delta\Gamma=\delta_0$」的分布语言定义了基本解，这里则给出了它的 Fourier 侧对应 $\widehat\Gamma(\xi)=|\xi|^{-2}$。

**补充（热半群的 Fourier 侧）**　由 21.4.2，$-\Delta$ 在 Fourier 侧就是乘 $|\xi|^2$；于是 [[20-Duhamel原理与齐次化]] §20.5 里那个抽象半群 $e^{t\Delta}$，在 Fourier 侧就是**逐点乘 $e^{-|\xi|^2t}$**。这正是 [[22-Cauchy问题与热核]] 的全部内容。

---

## 订正

1. **Plancherel 证明里的辅助恒等式左端漏了共轭。** 讲稿写
$\int\mathcal F[f](\xi)h(\xi)d\xi=(2\pi)^n\int f(x)\overline{\mathcal F^{-1}[h]}(x)dx$，
左端应为 $\int\mathcal F[f]\overline{h}\,d\xi$。**判据**：若按原式（左端无共轭）取 $h=\overline{\mathcal F[g]}$，右端给出 $(2\pi)^n\int fg$，而定理要的是 $(2\pi)^n\int f\bar g$ —— 差一个共轭，只有在 $g$ 实值时才对。讲稿接着的计算里第二行是 $\int h(\xi)e^{-\sqrt{-1}\langle x,\xi\rangle}d\xi$、第三行却变成 $\int\overline{h(\xi)}e^{\sqrt{-1}\langle\xi,x\rangle}d\xi$，共轭是**中途冒出来**的，说明原意就是 $\bar h$。本篇 §21.5 已按 $\bar h$ 重写。
2. **§5 里符号 $u_\xi$ 从未定义。** 讲稿写「$\varphi(x)=\frac{1}{(2\pi)^n}\int\widehat\varphi(\xi)u_\xi\,d\xi$」，但 $u_\xi$ 前文只以 $e^{\sqrt{-1}\langle\xi,x\rangle}$ 的形式出现过（且是在讲 $-\Delta$ 的特征函数时），没有赋名。应为 $u_\xi(x)=e^{\sqrt{-1}\langle\xi,x\rangle}$，且这里必须写 $u_\xi(x)$ 而不是 $u_\xi$（它是 $x$ 的函数）。
3. **「此时特征值 $\{|\xi|^2\}=\mathbb R_+$」有两处不准确。** 其一，$|\xi|^2$ 取值于 $[0,+\infty)$，包含 $0$（$\xi=0$，对应常数函数），而 $\mathbb R_+$ 常指 $(0,+\infty)$；其二，连续情形下 $|\xi|^2$ 不是「特征值」（$u_\xi\notin L^2$），而是**谱点**（连续谱）。正确说法：$-\Delta$ 在 $\mathbb R^n$ 上的谱是 $[0,+\infty)$，且全部是连续谱、没有特征值。
4. **定理 21.4.1 的标题说「同构」，证明只证了「映进」。** 见 §21.4.4 的补充：另一半由 $\mathcal F^{-1}[g](x)=(2\pi)^{-n}\mathcal F[g](-x)$ 与 $\mathcal S$ 对反射封闭得到，讲稿未写。
5. **卷积定理的证明跳过了「为什么可以 Fubini」。** $f_1,f_2\in\mathcal S$ 使 $|f_1(x-y)f_2(y)|$ 关于 $(x,y)$ 可积（两个速降函数之积在 $\mathbb R^{2n}$ 上可积），这是交换积分次序的依据；讲稿只写「其中第二个等号由 Fubini 定理得到」，没给可积性。本篇已补一句。
6. **反演公式的证明里，标着「$f\in\mathcal S$」的那次交换积分次序其实不合法。** 讲稿第九周写

$\mathcal F^{-1}[\mathcal F[f]](y)=\frac{1}{(2\pi)^n}\int_{\mathbb R^n}\int_{\mathbb R^n}f(x)e^{-\sqrt{-1}\langle x,\xi\rangle+\sqrt{-1}\langle y,\xi\rangle}dx\,d\xi\stackrel{f\in\mathcal S}{=}\frac{1}{(2\pi)^n}\int_{\mathbb R^n}\left[\int_{\mathbb R^n}e^{\sqrt{-1}\langle y-x,\xi\rangle}d\xi\right]f(x)\,dx$，

但 Fubini 要求 $\iint|f(x)|\,dx\,d\xi<+\infty$，而 $\int_{\mathbb R^n}\left|e^{\sqrt{-1}\langle y-x,\xi\rangle}\right|d\xi=+\infty$ —— **内层积分根本不收敛**，与 $f$ 是不是速降的毫无关系（$f$ 只控制 $x$ 方向）。这一步只能是形式推导，**真正的证明从下面的 $I_\varepsilon$ 开始**，那个 $e^{-\frac\varepsilon2|\xi|^2}$ 就是用来让内层积分收敛的。讲稿把它写成一个带依据的等号，读者会以为速降性足够，从而错过扰动法存在的理由。本篇 §21.3.1 已把这一步明确标为「不合法，只是形式」，并把证明改写成「一个量、两种算法、各取极限」。
7. **扰动法中间一步的变量混用。** 讲稿写

$I_\varepsilon(y)=\frac{1}{\pi^{n/2}}\int_{\mathbb R^n}e^{-\frac{|x|^2}{2\varepsilon}}f\left(y-\sqrt{2\varepsilon}\frac{x}{\sqrt{2\varepsilon}}\right)d\left(\frac{x}{\sqrt{2\varepsilon}}\right)$，

把旧变量 $x$ 与新变量 $\zeta=x/\sqrt{2\varepsilon}$ 混在同一行：指数上用的是旧 $x$（应换成 $\zeta$ 并约掉 $\varepsilon$），$f$ 里却写成 $\sqrt{2\varepsilon}\cdot\frac{x}{\sqrt{2\varepsilon}}=x$（等于没换元），而微分已经换成了 $d\zeta$。结果 $e^{-\frac{|x|^2}{2\varepsilon}}$ 的 $\varepsilon$ 没有消掉，与末行 $e^{-|\xi|^2}$ 对不上。正确写法只有一个新变量：令 $x=\sqrt{2\varepsilon}\zeta$，得 $\frac{1}{\pi^{n/2}}\int e^{-|\zeta|^2}f(y-\sqrt{2\varepsilon}\zeta)\,d\zeta$（本篇 (21.7)）。

---

系列：[[20-Duhamel原理与齐次化]] · [[22-Cauchy问题与热核]] · [[15-Sturm-Liouville定理与一维特征值问题]]
