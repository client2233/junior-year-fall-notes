# 20 多项式逼近·三角逼近与离散 Fourier 变换

> 接 [[19-多项式逼近·离散最小二乘]]。多项式逼近有一个结构性缺陷：它对**周期函数**不合适。用多项式去逼近 $\sin x$，区间一长次数就爆掉（多项式必然在远处发散，而周期函数不会）。
>
> 正确的基是三角函数。把 [[17-多项式逼近·内积空间的最佳逼近]] 的框架换一组正交基，得到的就是 Fourier 展开；再离散化，就是 DFT —— 整章最后一个结果，也是最有用的一个。

## 3.7 周期函数类的最佳平方逼近与离散 Fourier 变换

### 3.7.1 三角多项式空间

设 $L^{2}([0,2\pi])$ 为周期 $2\pi$ 的实值平方可积函数空间，内积与范数为

$$
\langle f,g\rangle=\int_{0}^{2\pi}f(x)g(x)\,\mathrm dx,\qquad \lVert f\rVert=\sqrt{\langle f,f\rangle}.
$$

取三角多项式空间

$$
X=\operatorname{span}\{\varphi_{0},\varphi_{1},\dots,\varphi_{2n}\},
$$

其中

$$
\varphi_{0}(x)=\frac12,\qquad \varphi_{2k}(x)=\cos kx,\qquad \varphi_{2k-1}(x)=\sin kx,\qquad k=1,2,\dots,n.
$$

$\dim X=2n+1$，$\{\varphi_{j}\}_{j=0}^{2n}$ 是 $X$ 的一组**正交基**。

**验正交性.** 全部靠三个基本积分（$k,m\ge1$）：

$$
\int_{0}^{2\pi}\cos kx\,\mathrm dx=0,\qquad \int_{0}^{2\pi}\sin kx\,\mathrm dx=0,
$$

$$
\int_{0}^{2\pi}\cos kx\cos mx\,\mathrm dx=\pi\delta_{km},\qquad \int_{0}^{2\pi}\sin kx\sin mx\,\mathrm dx=\pi\delta_{km},\qquad \int_{0}^{2\pi}\cos kx\sin mx\,\mathrm dx=0 .
$$

因此

$$
\langle\varphi_{0},\varphi_{0}\rangle=\int_{0}^{2\pi}\frac14\,\mathrm dx=\frac{\pi}{2},\qquad
\langle\varphi_{j},\varphi_{j}\rangle=\pi\ (j\ge1),\qquad
\langle\varphi_{i},\varphi_{j}\rangle=0\ (i\ne j).
$$

> **要点（$\varphi_{0}=\frac12$ 这个奇怪的选择）** 为什么不取 $\varphi_{0}=1$？因为那样 $\langle\varphi_{0},\varphi_{0}\rangle=2\pi$，而其余基的范数平方是 $\pi$ —— 差 2 倍，后面每个公式都要带一个例外项。
>
> 取 $\varphi_{0}=\frac12$ 后 $\langle\varphi_{0},\varphi_{0}\rangle=\frac{\pi}{2}$，**仍然**不是 $\pi$。真正的效果在系数公式上：见 3.7.2，$a_{0}$ 与 $a_{k}$ 能用同一个公式 $\frac{1}{\pi}\int f\cos kx\,\mathrm dx$（$k=0$ 时 $\cos0=1$），代价是展开式里 $a_{0}$ 那一项要写成 $\frac{a_{0}}{2}$。**这是标准 Fourier 级数记号「$a_{0}/2+\sum$」的真正来源** —— 不是习惯，是为了让系数公式统一。

### 3.7.2 最佳平方逼近 = Fourier 级数部分和

由最佳平方逼近元的存在唯一性（[[17-多项式逼近·内积空间的最佳逼近]] 定理 3.8），对任意周期 $2\pi$ 的 $f$，存在唯一的 $f_{n}\in X$ 使

$$
\lVert f_{n}-f\rVert=\inf_{g\in X}\lVert g-f\rVert.
$$

由特征定理，$\langle f-f_{n},\varphi_{j}\rangle=0$，$j=0,\dots,2n$，逐项写出：

$$
f_{n}(x)=\frac{a_{0}}{2}+\sum_{k=1}^{n}\bigl(a_{k}\cos kx+b_{k}\sin kx\bigr),
$$

其中

$$
a_{k}=\frac{\langle f,\varphi_{2k}\rangle}{\langle\varphi_{2k},\varphi_{2k}\rangle}=\frac{1}{\pi}\int_{0}^{2\pi}f(x)\cos kx\,\mathrm dx,\qquad k=0,1,\dots,n,
$$

$$
b_{k}=\frac{\langle f,\varphi_{2k-1}\rangle}{\langle\varphi_{2k-1},\varphi_{2k-1}\rangle}=\frac{1}{\pi}\int_{0}^{2\pi}f(x)\sin kx\,\mathrm dx,\qquad k=1,2,\dots,n.
$$

> **要点（两个「1/π」为什么统一）** $k\ge1$ 时 $\langle\varphi_{2k},\varphi_{2k}\rangle=\pi$，所以系数是 $\frac{1}{\pi}\int$。$k=0$ 时 $\langle\varphi_{0},\varphi_{0}\rangle=\frac{\pi}{2}$，系数应该是 $\frac{2}{\pi}\int\frac12 f=\frac{1}{\pi}\int f$ —— **也是 $\frac{1}{\pi}\int$**。所以 $a_{k}$ 的公式对 $k=0$ 自动成立，不需要单独定义。
>
> 但展开式里不能写 $a_{0}\varphi_{0}=a_{0}\cdot\frac12$ 之外的东西 —— 它就是 $\frac{a_{0}}{2}$。所以「$a_{0}/2+\cdots$」不是笔误，是 $\varphi_{0}=\frac12$ 的必然结果。

> **要点（$f_{n}$ 就是 Fourier 级数的部分和）** 线性代数的话说：$f_{n}$ 是 $f$ 在 $X$ 上的**正交投影**。分析学的话说：$f_{n}$ 是 $f$ 的 Fourier 级数前 $n$ 项和。$a_{k},b_{k}$ 就是 Fourier 系数。
>
> 两种说法完全等价，但**线性代数的说法更好用** —— 因为它顺便告诉你：$f_{n}$ 是 $X$ 中**唯一**的最佳 $L^{2}$ 逼近（不需要讨论 「Fourier 级数是否收敛」这个问题）。

### 3.7.3 平方收敛与一致收敛

由 Fourier 级数的收敛性定理：若 $f$ 在 $[0,2\pi]$ 上平方可积，则

$$
\lim_{n\to\infty}\lVert f_{n}-f\rVert=0\qquad\text{（}L^{2}\text{ 收敛）}.
$$

> **要点（$L^{2}$ 收敛 $\not\Rightarrow$ 一致收敛）** 这是本篇最重要的警告。$L^{2}$ 收敛只说「误差的能量趋于零」，完全允许误差在个别点附近尖得很高（能量小但峰值大）。
>
> 经典反例是方波：$f(x)=\operatorname{sign}(\sin x)$ 的 Fourier 部分和在跳跃点附近始终有约 $9\%$ 的过冲（**Gibbs 现象**），$n\to\infty$ 也不消失。它的 $L^{2}$ 误差确实趋于零，但 $\lVert f_{n}-f\rVert_{\infty}$ 不趋于零。
>
> 一致收敛需要 $f$ 更好的性质（例如 $f$ 连续且导数平方可积，或 $f$ 属于 Hölder 类）。对照 [[18-多项式逼近·最佳平方逼近与正交多项式]] 3.5.8 的 Chebyshev 情形：那里 $f\in C^{1}$ 就能保证一致收敛，因为权函数在端点发散起到了「压制边界振荡」的作用。Fourier 情形没有这个机制，跳跃点就是跳跃点。

### 3.7.4 离散情形：复三角多项式拟合

**复内积.** 考虑 $\mathbb C^{n}$，内积与范数为

$$
\langle\beta,\gamma\rangle=\sum_{k=1}^{n}\beta_{k}\overline{\gamma_{k}},\qquad \lVert\gamma\rVert_{\mathbb C^{n}}=\sqrt{\langle\gamma,\gamma\rangle},
$$

其中 $\overline{\gamma_{k}}$ 是复共轭。

> **要点（复内积对第二变量共轭线性）** 这个定义满足 $\langle\beta,\gamma\rangle=\overline{\langle\gamma,\beta\rangle}$，且对**第一**变量线性、对**第二**变量共轭线性。实内积的对称性 $\langle x,y\rangle=\langle y,x\rangle$ 在这里升级成共轭对称。
>
> 这不是形式主义：3.7.5 的正交性计算里，$e^{i(j-\ell)x_{k}}$ 这个漂亮的结果完全依赖 $\overline{e^{i\ell x}}=e^{-i\ell x}$。如果误用实内积的定义（不取共轭），$\langle X_{j},X_{\ell}\rangle=\sum e^{i(j+\ell)x_{k}}$，正交性立刻消失。

**空间与基.** 记 $T_{m}$ 为所有次数不超过 $m-1$ 的复系数三角多项式构成的空间，则 $T_{m}$ 有一组基

$$
\varphi_{j}(x)=\mathrm{e}^{ijx},\qquad j=0,1,\dots,m-1,\qquad \dim T_{m}=m .
$$

**数据与问题.** 设数据为

$$
(x_{j},y_{j}),\qquad x_{j}=\frac{2\pi(j-1)}{n},\quad j=1,2,\dots,n,
$$

即 $\{x_{j}\}$ 是 $[0,2\pi)$ 上的 $n$ 个**等距**点。要在 $T_{m}$ 中（$m\le n$）找

$$
s^{*}(x)=\sum_{j=0}^{m-1}\alpha_{j}^{*}\mathrm{e}^{ijx}
$$

使

$$
\sum_{j=1}^{n}\bigl|s^{*}(x_{j})-y_{j}\bigr|^{2}=\inf_{s\in T_{m}}\sum_{j=1}^{n}\bigl|s(x_{j})-y_{j}\bigr|^{2}.
$$

称 $s^{*}$ 为数据 $\{(x_{j},y_{j})\}$ 的**最小二乘 $m-1$ 次三角多项式**。

### 3.7.5 等距点上的正交性

把问题转成 $\mathbb C^{n}$ 中的最佳逼近。令

$$
X_{j}=\bigl(\varphi_{j}(x_{1}),\dots,\varphi_{j}(x_{n})\bigr)^{T}=\bigl(\mathrm{e}^{i(j-1)x_{1}},\dots,\mathrm{e}^{i(j-1)x_{n}}\bigr)^{T},\qquad j=1,\dots,m .
$$

利用 $x_{k}=\frac{2\pi(k-1)}{n}$ 与**等比数列求和**：

$$
\langle X_{j},X_{\ell}\rangle_{\mathbb C^{n}}=\sum_{k=1}^{n}\mathrm{e}^{i(j-\ell)x_{k}}=\sum_{k=0}^{n-1}\mathrm{e}^{i(j-\ell)\frac{2\pi k}{n}}=
\begin{cases}
0, & j\ne\ell,\\[2pt]
n, & j=\ell,
\end{cases}
$$

（当 $j\ne\ell$ 时公比 $\omega=\mathrm{e}^{i(j-\ell)2\pi/n}\ne1$ 且 $\omega^{n}=1$，故和为 $\frac{1-\omega^{n}}{1-\omega}=0$。）

所以 $\{X_{j}\}_{j=1}^{m}$ 是 $\mathbb C^{n}$ 中的**正交**向量组 —— 而且**与 $n$ 无关、与区间无关，只要等距就成立**。

> **要点（这条正交性是整个 DFT 的根基）** 换成非等距点，$\langle X_{j},X_{\ell}\rangle$ 一般不为零，最小二乘就要解 $m\times m$ 的方程组。**等距采样把 Gram 矩阵变成了 $nI$**，于是系数一个除法就出来。
>
> 这条性质的离散版本，与 3.7.1 里 $\int\cos kx\cos mx=\pi\delta_{km}$ 的连续版本是同一件事的不同采样。**「等距」是保证它们一致的唯一条件** —— 这就是为什么采样定理要求等间隔。

### 3.7.6 系数公式与 DFT

令 $\Lambda=(\alpha_{1},\dots,\alpha_{m})^{T}$、$Y=(y_{1},\dots,y_{n})^{T}$、矩阵 $P=(X_{1},\dots,X_{m})\in\mathbb C^{n\times m}$。问题化为

$$
\lVert P\Lambda^{*}-Y\rVert_{\mathbb C^{n}}=\inf_{\Lambda\in\mathbb C^{m}}\lVert P\Lambda-Y\rVert_{\mathbb C^{n}}.
$$

由 $\{X_{j}\}$ 的正交性，

$$
\lambda_{j}^{*}=\frac{1}{n}\langle Y,X_{j}\rangle_{\mathbb C^{n}}=\frac{1}{n}\sum_{k=0}^{n-1}y_{k+1}\,\mathrm{e}^{-i(j-1)2\pi k/n},\qquad j=1,2,\dots,m. \tag{3.8}
$$

**当 $m=n$ 时退化为插值.** 此时 $s^{*}(x)=\sum_{j=0}^{n-1}\alpha_{j}^{*}\mathrm{e}^{ijx}$ 的自由度等于数据个数，

$$
\bigl(s^{*}(x_{1}),\dots,s^{*}(x_{n})\bigr)^{T}=Y,\qquad\text{即}\quad s^{*}(x_{j})=y_{j},\ j=1,\dots,n .
$$

也就是说

$$
y_{j}=\sum_{k=0}^{n-1}\lambda_{k+1}^{*}\,\mathrm{e}^{i2\pi k(j-1)/n},\qquad j=1,2,\dots,n, \tag{3.9}
$$

其中

$$
\lambda_{k}^{*}=\frac{1}{n}\sum_{j=0}^{n-1}y_{j+1}\,\mathrm{e}^{-i(k-1)2\pi j/n},\qquad k=1,2,\dots,n. \tag{3.10}
$$

> **要点（(3.9)(3.10) 就是 DFT 与 IDFT）** 把下标平移一下（$j\to j+1$，$k\to k+1$，(3.10) 里的 $(k-1)$ 就变成 $k$），得到标准形式：

$$
\boxed{\ \hat a_{j}=\sum_{\ell=0}^{N-1}a_{\ell}\,\mathrm{e}^{-2\pi ij\ell/N},\qquad a_{j}=\frac{1}{N}\sum_{\ell=0}^{N-1}\hat a_{\ell}\,\mathrm{e}^{2\pi ij\ell/N}\ }
$$

> 前者是**离散 Fourier 变换**（DFT），后者是**逆变换**（IDFT）。
>
> 两者的差别只有两点：**指数符号相反**、**IDFT 多一个 $\frac{1}{N}$**。这个 $1/N$ 的来源不是约定，而是 $\lVert X_{j}\rVert^{2}=n$ —— 正交但不是标准正交，所以要除以范数平方。

**对应关系.** 在最小二乘问题中取 $N=n$、$\hat a_{j}=y_{j+1}$、$a_{j}=\lambda_{j+1}^{*}$，就看出**离散最佳平方逼近问题与 DFT/IDFT 完全是一回事**。

> **要点（一句话总结全章）** 从 $\S3.1$ 到这里，所有逼近问题都是「在有限维子空间中求最近点」。Fourier 变换只是把子空间取成三角多项式、把点取成等距时的特例。**DFT 不是新东西，是正交投影在等距采样下的闭式解。**

### 3.7.7 数值例子

$n=6$，$y=(1,\,2,\,0,\,-1,\,3,\,0.5)$：

| $j$ | $\hat a_{j}$（DFT） | $a_{j}=\hat a_{j}/6$ |
|---|---|---|
| 0 | $5.5$ | $0.916667$ |
| 1 | $1.75+1.29904i$ | $0.291667+0.216506i$ |
| 2 | $-2.75-3.89711i$ | $-0.458333-0.649519i$ |
| 3 | $2.5$ | $0.416667$ |
| 4 | $-2.75+3.89711i$ | $-0.458333+0.649519i$ |
| 5 | $1.75-1.29904i$ | $0.291667-0.216506i$ |

三个自检：

1. **共轭对称**：$\hat a_{6-j}=\overline{\hat a_{j}}$（因 $y$ 实）。这正是「$n$ 个实数只有 $n$ 个实自由度」的体现。
2. **IDFT 复原**：$\frac{1}{6}\sum_{j}\hat a_{j}\mathrm{e}^{2\pi ij\ell/6}$ 逐位回到 $(1,2,0,-1,3,0.5)$。
3. **Parseval**：$\sum_{j=0}^{5}|y_{j}|^{2}=15.25$，$\frac{1}{6}\sum_{j=0}^{5}|\hat a_{j}|^{2}=\frac{91.5}{6}=15.25$ ✓

> **要点（Parseval 是投影的长度守恒）** 这条等式就是 [[17-多项式逼近·内积空间的最佳逼近]] 里的 Parseval 等式在有限维、完备正交系（$m=n$ 时 $\{X_{j}\}$ 张满 $\mathbb C^{n}$）下的形态。因为 $m=n$ 是插值，误差为零，Bessel 不等式取等。

> **补充（FFT 的位置）** 直接按定义算 DFT 是 $O(n^{2})$ 次复数乘加。Cooley–Tukey 的 FFT 利用 $\mathrm{e}^{-2\pi i/N}$ 的周期性把计算量降到 $O(n\log n)$（$n=2^{20}$ 时约差 $5\times10^{4}$ 倍）。
>
> 关键点：**FFT 算的仍然是定义 (3.7.6) 里的那个 DFT**，一点精度都没牺牲（理想算术下逐位相同）。它的作用纯粹是让「$n$ 很大」从不可能变成日常 —— 音频、图像、通信里的 Fourier 方法之所以实用，全靠这一个算法。
>
> 顺带一句：DFT 矩阵 $F=(\mathrm{e}^{-2\pi ij\ell/N})_{j,\ell}$ 满足 $\frac{1}{\sqrt N}F$ 是**酉矩阵**（$F^{*}F=NI$）。这就是为什么 IDFT 只需除以 $N$、不需要别的校正 —— 酉性保证了变换可逆且长度守恒。

### 3.7.8 本节小结

- 三角多项式空间以 $\{\frac12,\cos kx,\sin kx\}$ 为正交基；$\varphi_{0}=\frac12$ 的选择使 $a_{k}$ 的公式对 $k=0$ 统一，代价是展开式的 $\frac{a_{0}}{2}$。
- 最佳平方逼近 = 正交投影 = Fourier 级数部分和。存在唯一性来自有限维内积空间理论，不需要另外讨论收敛。
- $L^{2}$ 收敛 $\not\Rightarrow$ 一致收敛（Gibbs 现象）。一致收敛要 $f$ 更光滑。
- 复内积必须对第二变量取共轭，否则等距点上的正交性不成立。
- 等距采样使 $\langle X_{j},X_{\ell}\rangle=n\delta_{j\ell}$，Gram 矩阵化成 $nI$，系数闭式求出 —— 这就是 DFT。
- $m=n$ 时最小二乘退化为插值，反解就是 IDFT。DFT/IDFT 只差指数符号与 $1/N$。

## 本节索引

- 内积空间、特征定理、正交基、Bessel/Parseval ← 本篇的全部框架，逐条平移 ← [[17-多项式逼近·内积空间的最佳逼近]]
- 正交基与 Gram–Schmidt ← 三角基的正交性验证方式 ← [[18-多项式逼近·最佳平方逼近与正交多项式]]
- 离散内积、最小二乘、等距数据 → 3.7.4 的直接前置 ← [[19-多项式逼近·离散最小二乘]]
- Chebyshev 展开（$x=\cos\theta$ 换元） → 与 Fourier 展开的对照：两者都是把函数展成三角函数的组合，只是一个在 $\theta$ 上、一个在 $x$ 上 ← [[18-多项式逼近·最佳平方逼近与正交多项式]]
- $L^{2}$ 收敛与一致收敛的差别 → $C[a,b]$ 上的一致范数作对照 ← [[13-多项式逼近·赋范空间与最佳逼近元的存在性]]
- Parseval 等式（泛函分析版） ← [[11-Hilbert空间·正交基]]

---

> **本章完。** 第三章从「逼近」这个提法出发，走完了一致逼近（交错定理 + Remes）与平方逼近（正交投影 + 正交多项式 + 最小二乘 + DFT）两条线。下一章进入数值积分。

---

系列：[[19-多项式逼近·离散最小二乘]]
