# 16 特征值的变分刻画与 Poincaré 不等式

> 系列：[[15-Sturm-Liouville定理与一维特征值问题]] · [[14-谱理论·鼓面的数学建模与Laplace特征值]] · [[17-矩形与圆盘上的分离变量·Bessel方程]]

[[14-谱理论·鼓面的数学建模与Laplace特征值]] 证了谱的存在与离散性，但没说 $\lambda_1$ **等于什么**。本篇给三种刻画：

1. **能量恒等式**：把三类边值统一成「$\lambda$ = 某个商」的形式，顺手得 Robin 的 $\lambda>0$；
2. **Rayleigh 商**：$\lambda_1=\inf\int|\nabla u|^2/\int u^2$ —— 不用解方程就能估计第一特征值；
3. **Poincaré 不等式**：前面那个商的倒数就是最佳常数，并且有一个只依赖 $\operatorname{diam}\Omega$ 的下界。

这三条是「不知道 $\Omega$ 长什么样时也能说点话」的全部本钱，后面的热方程衰减率、数值方法的误差估计都要用它。

---

## 16.1 Robin 条件下的非负性

Dirichlet 与 Neumann 的 $\lambda\ge0$ 已在 [[14-谱理论·鼓面的数学建模与Laplace特征值]] §14.3 证完，这里补 Robin。设

$$
\begin{cases}-\Delta u=\lambda u, & x\in\Omega,\\[4pt] \dfrac{\partial u}{\partial\nu}+\sigma u=0, & x\in\partial\Omega,\end{cases}\qquad\sigma>0 .
$$

用 Green 第一恒等式把 $\int u(-\Delta u)$ 拆开：

$$
\int_\Omega u\cdot(-\Delta u)\,dx=-\int_{\partial\Omega}\frac{\partial u}{\partial\nu}\cdot u\,d\sigma_x+\int_\Omega|\nabla u|^2\,dx .
$$

边界上 $\dfrac{\partial u}{\partial\nu}=-\sigma u$，故第一项为 $\sigma\int_{\partial\Omega}u^2\,d\sigma_x$：

$$
\int_\Omega u\cdot(-\Delta u)\,dx=\sigma\int_{\partial\Omega}u^2\,d\sigma_x+\int_\Omega|\nabla u|^2\,dx .
$$

由方程，左端又等于 $\lambda\int_\Omega u^2\,dx$。于是

$$
\lambda=\frac{\sigma\displaystyle\int_{\partial\Omega}u^2\,d\sigma_x+\int_\Omega|\nabla u|^2\,dx}{\displaystyle\int_\Omega u^2\,dx}>0 .
$$

**不仅 $\ge0$，而且严格大于 $0$**：分子含 $\sigma\int_{\partial\Omega}u^2$ 这一项，若它为零则 $u|_{\partial\Omega}=0$（$u$ 连续），于是 $u$ 是 Dirichlet 特征函数、$\lambda>0$；否则该项严格为正。

> **要点（为什么 Robin 夹在 Dirichlet 与 Neumann 之间）**
> 上式右端是个商，分子 = $\int|\nabla u|^2+\sigma\int_{\partial\Omega}u^2$。对固定的 $u$，这个分子随 $\sigma$ 单调递增，所以 $\sigma\mapsto\lambda_1(\sigma)$ 单调不减；$\sigma\to0^+$ 时退化为 Neumann（$\lambda_1\to0$），$\sigma\to+\infty$ 时边界项强迫 $u|_{\partial\Omega}\to0$、退化为 Dirichlet。后面 [[17-矩形与圆盘上的分离变量·Bessel方程]] 与 [[18-热方程·建模、边值条件与分离变量法]] 里出现的「$\lambda_1$ 决定衰减率」，说的都是这个商的下确界。

---

## 16.2 Rayleigh 商：$\lambda_1$ 的变分刻画

**定理 16.2.1（Rayleigh 商）**　在 Dirichlet 边值条件下，

$$
\lambda_1=\inf_{\substack{u\in C^\infty(\overline\Omega)\\ u|_{\partial\Omega}=0,\ u\not\equiv0}}\frac{\displaystyle\int_\Omega|\nabla u|^2\,dx}{\displaystyle\int_\Omega u^2\,dx}.
$$

**证明**　记右端的分式为 $R(u)$（**Rayleigh 商**）。

**下界（$R(u)\ge\lambda_1$）。** 任取 $u\in C^\infty(\overline\Omega)$、$u|_{\partial\Omega}=0$。由谱定理，$\{u_k\}$ 是 $L^2(\Omega)$ 的完备标准正交基，故

$$
u=\sum_{i=1}^{\infty}\alpha_iu_i,\qquad \alpha_i=\int_\Omega uu_i\,dx,\qquad \sum_{i=1}^{\infty}\alpha_i^2=\lVert u\rVert_{L^2}^2 .
$$

先用 Green 第一恒等式把梯度积分换成 $\int u(-\Delta u)$（边界项因 $u|_{\partial\Omega}=0$ 消失）：

$$
\int_\Omega|\nabla u|^2\,dx=\int_\Omega u\cdot(-\Delta u)\,dx .
$$

再把 $-\Delta u=\sum_i\alpha_i\lambda_iu_i$ 与 $u=\sum_j\alpha_ju_j$ 代入，用正交性：

$$
\int_\Omega|\nabla u|^2\,dx=\sum_{i=1}^{\infty}\alpha_i\lambda_i\int_\Omega u\,u_i\,dx=\sum_{i=1}^{\infty}\lambda_i\alpha_i^2\ge\lambda_1\sum_{i=1}^{\infty}\alpha_i^2=\lambda_1\lVert u\rVert_{L^2}^2 .
$$

故 $R(u)\ge\lambda_1$，即 $\lambda_1\le\inf R(u)$。

**上界（$\inf R(u)\le\lambda_1$）。** 取 $u=u_1$（第一特征函数），此时 $\int|\nabla u_1|^2=\lambda_1\int u_1^2$，故 $R(u_1)=\lambda_1$。既然下确界不超过任何一个取值，$\inf R\le\lambda_1$。

两向合并即得等号。$\square$

> **要点（这一步为什么重要）**　定理把「解 PDE」换成了「在一个函数空间上取下确界」。这有三重好处：
> 1. **估上界只要试函数**：任意取一个满足边值的 $u$ 代进 $R(u)$，就得到 $\lambda_1$ 的一个上界 —— 这是变分法、有限元、Rayleigh–Ritz 方法的全部起点；
> 2. **估下界要靠不等式**：比如 §16.3 的 Poincaré 常数，就是给 $R(u)$ 找一个只依赖 $\Omega$ 几何量的正下界；
> 3. **特征函数就是极小元**：$\lambda_1$ 的那个「$\inf$」在 $u=u_1$ 处取到，所以 $u_1$ 是 $R$ 的极小点，令一阶变分为零就回到特征值方程 —— 两者等价。

**补充（严格论证要换空间）**　上面的推导里有两处需要交代：(i) 无穷级数能否逐项求导 / 逐项积分；(ii) 为什么 $\inf$ 要在 $C^\infty(\overline\Omega)$ 而不是在 $W^{1,2}_0(\Omega)=H^1_0(\Omega)$ 中取。标准处理是：

- 把 $R$ 定义在 $H^1_0(\Omega)\setminus\{0\}$ 上，由 $H^1$ 的完备性 + 下确界 + 弱紧性（Rellich 紧嵌入）可证下确界在某个 $u_1\in H^1_0$ 处取到；
- 再对 $u_1$ 用椭圆正则性得 $u_1\in C^\infty(\Omega)$，边界光滑时甚至 $u_1\in C^\infty(\overline\Omega)$；
- 于是在 $C^\infty$ 上取下确界与在 $H^1_0$ 上取下确界相等（$C^\infty_c$ 在 $H^1_0$ 中稠密）。

课堂按形式推导走，是因为结论本身不依赖这些细节；但在**数值方法**里必须用 $H^1_0$ 口径（有限元空间不是 $C^\infty$），这一点在 [[17-矩形与圆盘上的分离变量·Bessel方程]] 之后的课程会显出来。

**补充（$\lambda_k$ 也可以这样刻画，超出本课）**　同样的手法给出第 $k$ 个特征值的 **Courant–Fischer 极小极大**公式：

$$
\lambda_k=\min_{\substack{E\subseteq H^1_0(\Omega)\\ \dim E=k}}\ \max_{u\in E\setminus\{0\}}R(u),
$$

$k=1$ 就是定理 16.2.1。它说明「$\lambda_k$ 大」等于「任何 $k$ 维子空间里都有一个方向使 $R$ 很大」—— 这是 Weyl 渐近律 $\lambda_k\sim C_nk^{2/n}$ 的出发点。

---

## 16.3 Poincaré 不等式与它的最佳常数

**定理 16.3.1（Poincaré 不等式）**　存在常数 $C>0$，使对一切 $f\in C^1(\overline\Omega)$、$f|_{\partial\Omega}=0$，

$$
\int_\Omega|f|^2\,dx\le C\int_\Omega|\nabla f|^2\,dx .
$$

其**最佳**（最小）常数为

$$
C=\frac{1}{\lambda_1}.
$$

**证明**　**（$C=1/\lambda_1$ 可用）** 把定理 16.2.1 的 $R(u)\ge\lambda_1$ 两边取倒数（$\int|\nabla u|^2>0$，否则 $u$ 为常数、由边值条件为零）：

$$
\int_\Omega|f|^2\,dx\le\frac{1}{\lambda_1}\int_\Omega|\nabla f|^2\,dx .
$$

**（不能更小）** 取 $f=u_1$：由 $\int|\nabla u_1|^2=\lambda_1\int u_1^2$ 得 $\int u_1^2=\dfrac{1}{\lambda_1}\int|\nabla u_1|^2$，等号成立。故任何 $C<1/\lambda_1$ 都会在 $u_1$ 上失败。$\square$

> **要点**　「最佳常数 = $1/\lambda_1$」这句话的含义比它看起来强：它说**第一特征值就是 Poincaré 不等式的倒数**。所以「能不能估出 Poincaré 常数」和「能不能估出 $\lambda_1$」是同一个问题，而后者有大量的几何估计（下面的定理 16.3.2 是最粗的一个）。

**定理 16.3.2（只用直径的下界）**

$$
\lambda_1\ge\frac{1}{\left(\operatorname{diam}(\Omega)\right)^2},\qquad \operatorname{diam}(\Omega):=\sup_{x,y\in\Omega}|x-y| .
$$

**证明**　取 $e_1=(1,0,\dots,0)$，把 $u$（$u|_{\partial\Omega}=0$）零延拓到 $\mathbb R^n$。

**第一步：把 $u$ 写成沿 $e_1$ 的积分。** 对 $x\in\Omega$，

$$
u(x)=\int_{-\infty}^{0}\frac{d}{ds}u(x+se_1)\,ds=\int_{-\infty}^{0}\frac{\partial u}{\partial x_1}(x+se_1)\,ds ,
$$

因为 $s\to-\infty$ 时 $x+se_1$ 已在 $\Omega$ 外、$u$ 取零。又因 $x\in\Omega$ 时使 $x+se_1\in\Omega$ 的 $s\le0$ 只有一段、长度不超过 $\operatorname{diam}(\Omega)$，被积函数在该段之外为零，故可把积分截断到 $(-\operatorname{diam}(\Omega),0)$。

**第二步：Cauchy–Schwarz。** 对 $s$ 在长度不超过 $\operatorname{diam}$ 的区间上积分：

$$
|u(x)|\le\left(\int_{-\operatorname{diam}(\Omega)}^{0}\left|\nabla u\right|^2(x+se_1)\,ds\right)^{1/2}\left(\operatorname{diam}(\Omega)\right)^{1/2}.
$$

**第三步：平方并对 $x$ 积分。** 两边平方：

$$
u(x)^2\le\operatorname{diam}(\Omega)\cdot\int_{-\operatorname{diam}(\Omega)}^{0}|\nabla u|^2(x+se_1)\,ds .
$$

在 $\Omega$ 上对 $x$ 积分，再用 Fubini（交换 $x$ 与 $s$）：

$$
\int_\Omega u(x)^2\,dx\le\operatorname{diam}(\Omega)\int_{-\operatorname{diam}(\Omega)}^{0}\left[\int_\Omega|\nabla u|^2(x+se_1)\,dx\right]ds .
$$

内层对 $x\in\Omega$ 的积分可换成对 $x\in\mathbb R^n$ 的积分（$u$ 零延拓后 $\nabla u$ 在 $\Omega$ 外为零），每次换元 $x\mapsto x+se_1$ 不改变积分值：

$$
\int_\Omega|\nabla u|^2(x+se_1)\,dx\le\int_{\mathbb R^n}|\nabla u|^2=\int_\Omega|\nabla u|^2 .
$$

外层对 $s$ 的区间长度是 $\operatorname{diam}(\Omega)$，故

$$
\int_\Omega u^2\,dx\le\left(\operatorname{diam}(\Omega)\right)^2\int_\Omega|\nabla u|^2\,dx .
$$

由定理 16.3.1，最佳常数 $\dfrac{1}{\lambda_1}\le\left(\operatorname{diam}\Omega\right)^2$，即 $\lambda_1\ge\dfrac{1}{\left(\operatorname{diam}\Omega\right)^2}$。$\square$

> **要点（这个下界有多粗）**　一维情形 $\Omega=(0,L)$ 的精确值是 $\lambda_1=\pi^2/L^2$，而定理给的是 $1/L^2$ —— **差一个 $\pi^2$**。所以定理 16.3.2 只保证「量级对」，常数远远不是最优的。它的价值在于**完全不依赖 $\Omega$ 的形状**：不管多细长的香肠形区域，只要有界，$\lambda_1$ 就不会退化成 $0$。这个「形状无关」性质正是后面证明热方程解的唯一性、以及 $t\to\infty$ 时解指数衰减的关键。

---

## 订正

1. **Robin 情形的结论写成了 $\lambda\ge0$，应是 $\lambda>0$。** 讲稿 §16.1 由 $\lambda=\dfrac{\sigma\int_{\partial\Omega}u^2+\int_\Omega|\nabla u|^2}{\int_\Omega u^2}$ 推出 $\lambda\ge0$。但 $\sigma>0$ 时若分子为零，则 $u|_{\partial\Omega}=0$ 且 $\nabla u\equiv0$，$u\equiv0$，与特征函数的非零性矛盾 —— 所以分子恒正、$\lambda>0$。这与定理 15.2.1 的一维判别（$\lambda_1=0\iff$ Neumann）一致，写成 $\ge0$ 会漏掉「Robin 不可能有零特征值」这条信息。
2. **定理 16.2.1 的证明只做了形式推导，没交代级数的逐项求导 / 积分。** 讲稿直接写 $u=\sum\alpha_iu_i$ 后就对 $-\Delta u$ 展开。严格地说，需要 $u\in H^1_0$ 且 $\{u_i\}$ 是 $H^1$ 意义下的正交系（这由 $-\Delta u_i=\lambda_iu_i$ 与 Green 恒等式推出：$\int\nabla u_i\cdot\nabla u_j=\lambda_i\delta_{ij}$）。本篇在**补充**里给出了标准处理路线。
3. **「若 $u_1\in L^2(\Omega)$，由 $C^\infty$ 在 $W^{1,2}$ 中稠密性…」这段是防御性的。** 对 $\partial\Omega$ 光滑的有界区域，Dirichlet 特征函数由椭圆正则性其实属于 $C^\infty(\overline\Omega)$，所以「$u_1$ 只在 $L^2$ 里」这种情形不会发生。讲稿这段是为了在不假设边界正则性的前提下也能说通，读的时候不要误以为 $u_1$ 可能不光滑。
4. **定理 16.3.1 的「最佳」二字没给证明。** 讲稿只写「由特征值理论知 $C=1/\lambda_1$」。完整说法是两半：$C=1/\lambda_1$ 可用（由 Rayleigh 商）+ 不能更小（取 $f=u_1$ 取等）。缺了后半，「最佳」就没有依据。
5. **定理 16.3.2 的证明里 $\mathbb R^n$ 与 $\Omega$ 的积分混用。** 讲稿在
$\int_\Omega u^2\le\operatorname{diam}\int_{-\operatorname{diam}}^0\left[\int_\Omega|\nabla u|^2(x+se_1)dx\right]ds\le(\operatorname{diam})^2\int_{\mathbb R^n}|\nabla u|^2$
这几行里，内层先写 $\int_\Omega$、下一步写 $\int_{\mathbb R^n}$、最后一行又写回 $\int_\Omega$。三者相等是因为 $u$ 零延拓后 $\nabla u$ 在 $\Omega$ 外恒为零，但这一点必须点明，否则读者会以为换了积分区域。
6. **证明中「取 $e_1=(1,0,\dots,0)$」看起来像任意选定了一个方向，实际必须说明为什么一个方向就够。** 原因是对每个固定的 $x\in\Omega$，沿 $e_1$ 的射线段落在 $\Omega$ 内的长度都不超过 $\operatorname{diam}(\Omega)$ —— 而这个长度界与 $x$、与方向都无关，所以任取一个方向即可。讲稿把这句省了。

---

系列：[[15-Sturm-Liouville定理与一维特征值问题]] · [[14-谱理论·鼓面的数学建模与Laplace特征值]] · [[17-矩形与圆盘上的分离变量·Bessel方程]]
