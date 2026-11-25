# 17 矩形与圆盘上的分离变量·Bessel 方程

> 系列：[[16-特征值的变分刻画与Poincaré不等式]] · [[15-Sturm-Liouville定理与一维特征值问题]] · [[18-热方程·建模、边值条件与分离变量法]]

[[15-Sturm-Liouville定理与一维特征值问题]] 把一维问题算干净了。二维呢？**只有两种区域能算到底**：矩形与圆盘。它们能算，靠的都是「把二维算子拆成两个一维算子」—— 前者是对直角坐标拆，后者是对极坐标拆（拆出来的那一半不再是常系数方程，而是 **Bessel 方程**）。

讲稿这一节写得很短（两页），本篇把它铺开：补上分离变量的中间步骤、Bessel 方程的标准形式与换元、以及两个解的取舍理由。

---

## 17.1 矩形区域：谱的可加性

设 $\Omega=(0,a)\times(0,b)\subseteq\mathbb R^2$。取一维问题的解：$v=v(x)$ 满足

$$
-v''=\mu v,\quad x\in(0,a),\qquad w=w(y)\text{ 满足 } -w''=\rho w,\quad y\in(0,b),
$$

其中 $v,w$ 各带对应的边界条件（下面统一取 Dirichlet：$v(0)=v(a)=0$、$w(0)=w(b)=0$）。

令 $u(x,y):=v(x)w(y)$。直接算：

$$
-\Delta u(x,y)=-\bigl[v''(x)w(y)+v(x)w''(y)\bigr]=\mu v(x)w(y)+\rho v(x)w(y)=(\mu+\rho)\,u(x,y).
$$

而边界上（例如 $x=0$）$u(0,y)=v(0)w(y)=0$，四条边同理。所以

$$
(\mu+\rho,\ vw)\ \text{是}\ \Omega\ \text{上的特征对}.
$$

代入一维的显式谱（[[15-Sturm-Liouville定理与一维特征值问题]] §15.3）：

$$
\mu_m=\left(\frac{m\pi}{a}\right)^2,\quad m\ge1;\qquad \rho_n=\left(\frac{n\pi}{b}\right)^2,\quad n\ge1 .
$$

于是

$$
\lambda_{mn}=\pi^2\left(\frac{m^2}{a^2}+\frac{n^2}{b^2}\right),\qquad
u_{mn}(x,y)=\sin\frac{m\pi x}{a}\sin\frac{n\pi y}{b},\qquad m,n\ge1 .
$$

（归一化常数是 $2/\sqrt{ab}$，若要 $\lVert u_{mn}\rVert_{L^2}=1$ 需乘上它。）

> **要点（乘积给出的是全部特征函数吗）**　讲稿说「这是分离变量法的理论基础」，只说了一半：它证明了「乘积是特征函数」，**没证明「所有特征函数都是乘积」**。补齐靠一个纯测度论事实：若 $\{v_m\}$ 是 $L^2(0,a)$ 的完备标准正交基、$\{w_n\}$ 是 $L^2(0,b)$ 的，则 $\{v_m(x)w_n(y)\}_{m,n}$ 是 $L^2(\Omega)$ 的完备标准正交基（张量积积分的 Fubini 展开）。而由谱定理，$\Omega$ 上的特征函数也构成完备标准正交基；两族都是「由特征函数组成的完备正交基」，且特征值按重数计（$\dim E_\lambda$ 有限），故它们一致 —— **矩形上的特征函数恰好就是这些乘积**。没有这一步，「分离变量法」只是造出了一些解，不能说「全部」。

> **补充（简并与数论）**　正方形 $a=b$ 时，$\lambda_{mn}=\dfrac{\pi^2}{a^2}(m^2+n^2)$，于是 $(m,n)$ 与 $(n,m)$ 给出同一个特征值（$m\ne n$ 时重数 $\ge2$）。这种简并是**算术的**：$\lambda$ 的重数 = 方程 $m^2+n^2=N$ 的整数解个数，与「一个数有几种方式写成两个平方之和」的数论问题等价。圆形区域也有简并（见 §17.4），但来源不同 —— 那是旋转对称性造成的。

---

## 17.2 圆盘：极坐标下的分离变量

设 $B_1\subseteq\mathbb R^2$ 为单位圆盘，$(r,\theta)$ 为极坐标。由 [[12-习题课]] §12.1.2，极坐标下的 Laplace 算子是

$$
\Delta u=u_{rr}+\frac1ru_r+\frac{1}{r^2}u_{\theta\theta}.
$$

设 $u(x,y)=v(r)w(\theta)$，代入 $-\Delta u=\lambda u$：

$$
-\left[v''(r)w(\theta)+\frac1rv'(r)w(\theta)+\frac{1}{r^2}v(r)w''(\theta)\right]=\lambda v(r)w(\theta).
$$

两边除以 $vw$（$v,w$ 不恒为零），再乘 $r^2$：

$$
-\frac{r^2v''}{v}-\frac{rv'}{v}-\frac{w''}{w}=\lambda r^2 .
$$

移项，把只含 $r$ 的项与只含 $\theta$ 的项分到两边：

$$
-\frac{w''}{w}=\lambda r^2+\frac{r^2v''}{v}+\frac{rv'}{v}.
$$

左端只依赖 $\theta$、右端只依赖 $r$，故两边都等于同一个常数，记为 $k$。于是得到两个一维问题：

$$
w''+kw=0,
$$

$$
v''+\frac1rv'+\left(\lambda-\frac{k}{r^2}\right)v=0 .
$$

第一个方程的通解是 $w=A\cos\sqrt k\,\theta+B\sin\sqrt k\,\theta$。**几何上必须加上单值性条件**：$(r,\theta)$ 与 $(r,\theta+2\pi)$ 是同一点，故 $w$ 必须 $2\pi$-周期。这迫使 $k=m^2$、$m\in\mathbb Z_{\ge0}$，且

$$
w_m(\theta)=A_m\cos m\theta+B_m\sin m\theta .
$$

第二个方程就是所谓的 **Bessel 方程**（除以 $r^2$ 后的形式）：

$$
v''+\frac1rv'+\left(\lambda-\frac{m^2}{r^2}\right)v=0 .
$$

> **要点（$k=m^2$ 这一步不能省）**　讲稿直接写「$w$ 为三角函数」，中间跳了两步：先由「分离变量 ⇒ 常数为 $k$」，再由「单值性 ⇒ $k=m^2$ 且 $m$ 是整数」。若 $k<0$，$w$ 是双曲函数、不可能 $2\pi$-周期；若 $k>0$ 但不是完全平方，同样不周期。**周期性是圆盘问题独有的边值条件**（矩形没有），它是圆盘谱比矩形谱更「稀疏」的根源。

---

## 17.3 Bessel 方程与它的两个解

上式的系数含 $r$，不是常系数方程，但可以用**伸缩换元**化成标准形。令

$$
\tau=\sqrt\lambda\,r,\qquad v(r)=f(\tau).
$$

则 $v'(r)=\sqrt\lambda f'(\tau)$、$v''(r)=\lambda f''(\tau)$、$\dfrac1r=\dfrac{\sqrt\lambda}{\tau}$。代入：

$$
\lambda f''+\frac{\sqrt\lambda}{\tau}\cdot\sqrt\lambda f'+\left(\lambda-\frac{m^2\lambda}{\tau^2}\right)f=0 .
$$

三项都含因子 $\lambda$，约掉（$\lambda>0$）：

$$
f''+\frac1\tau f'+\left(1-\frac{n^2}{\tau^2}\right)f=0,\qquad n=m\in\mathbb Z_{\ge0} .\tag{17.1}
$$

这就是 **$n$ 阶 Bessel 方程**的标准形。它是二阶线性常微分方程，$\tau=0$ 是正则奇点，故有两个线性无关解：

$$
J_n(\tau)=\sum_{k=0}^{\infty}\frac{(-1)^k}{k!\,(k+n)!}\left(\frac{\tau}{2}\right)^{2k+n},
$$

$$
Y_n(\tau)=\frac{J_n(\tau)\cos n\pi-J_{-n}(\tau)}{\sin n\pi}\quad(n\notin\mathbb Z),\qquad
Y_n\ \text{在}\ n\in\mathbb Z\ \text{时取极限定义}.
$$

- **$J_n$**（第一类 Bessel 函数）在 $\tau=0$ 处有界：$J_0(0)=1$，$J_n(0)=0$（$n\ge1$），且 $J_n(\tau)\sim\dfrac{1}{n!}\left(\dfrac\tau2\right)^n$；
- **$Y_n$**（第二类，Neumann 函数）在 $\tau\to0^+$ 时发散：$Y_0(\tau)\sim\dfrac2\pi\ln\dfrac\tau2$，$Y_n(\tau)\sim-\dfrac{(n-1)!}{\pi}\left(\dfrac2\tau\right)^n$（$n\ge1$）。

**取舍**：圆盘内部（$r=0$）是定义域的内点，解必须有界，故

$$
v(r)=C\,J_m(\sqrt\lambda\,r),\qquad \text{舍去}\ Y_m .
$$

**定谱**：再代入 $r=1$ 的边值条件。

- **Dirichlet**（$v(1)=0$）：$J_m(\sqrt\lambda)=0$，即 $\sqrt\lambda$ 必须是 $J_m$ 的零点。记 $J_m$ 的正零点为

$$
0<j_{m,1}<j_{m,2}<\cdots\longrightarrow+\infty ,
$$

则

$$
\lambda_{m,k}=j_{m,k}^2,\qquad u_{m,k}(r,\theta)=J_m(j_{m,k}r)\bigl(A\cos m\theta+B\sin m\theta\bigr).
$$

- **Neumann**（$v'(1)=0$）：$J_m'(\sqrt\lambda)=0$，谱是 $J_m'$ 的零点平方。

> **补充（几个具体数值）**　$J_0$ 的前三个正零点为 $2.40483$、$5.52008$、$8.65373$；$J_1$ 的前两个为 $3.83171$、$7.01559$。于是单位圆盘 Dirichlet 问题的第一特征值 $\lambda_{0,1}=j_{0,1}^2\approx5.783$。**注意 $m\ge1$ 的每个 $\lambda_{m,k}$ 都是二重的**：$\cos m\theta$ 与 $\sin m\theta$ 张成二维特征空间 —— 这就是 §17.1 补充里说的「旋转对称导致的简并」，与矩形的数论简并来源不同。

---

## 17.4 两种区域对照

| | 矩形 $(0,a)\times(0,b)$ | 圆盘 $B_R$ |
|---|---|---|
| 坐标 | 直角 $(x,y)$ | 极 $(r,\theta)$ |
| 分离后 | 两个常系数方程 $-v''=\mu v$、$-w''=\rho w$ | $w''+m^2w=0$ 与 Bessel 方程 (17.1) |
| 特征函数 | 正弦/余弦的乘积 | $J_m(j_{m,k}r/R)\times$ 三角函数 |
| 谱 | $\pi^2(m^2/a^2+n^2/b^2)$，显式 | $j_{m,k}^2/R^2$，由零点定，**无闭式** |
| 简并来源 | 算术（$m^2+n^2$ 的表示数） | 旋转对称（$m\ge1$ 二重） |
| 节线 | 竖直/水平直线 | 同心圆（$J_m$ 的零点）+ 过圆心的直径（$\cos m\theta$ 的零点） |

> **要点**　**圆盘的谱算不出来，只能算到「零点」这一步** —— $j_{m,k}$ 是超越方程 $J_m(x)=0$ 的根，没有初等表达式。这与一维情形（$\lambda_n=(n\pi/L)^2$ 显式）形成鲜明对比。所以后面凡是需要圆盘谱的地方，要么用数值表，要么用 Weyl 渐近律 $\lambda_k\sim 4\pi k/|B_R|$ 这类定性结论。**矩形能算、圆盘只能定性，这个差别贯穿整门课。**

---

## 订正

1. **「该方程有两个线性无关的解：」后面是空的。** 讲稿在这一句之后直接跳到下一节（热方程），两个解根本没写出来。本篇已在 §17.3 补上 $J_n$ 与 $Y_n$，并写明「为什么舍掉 $Y_n$」——这一步是定谱的必要前提，漏掉它「$J_m(\sqrt\lambda)=0$」就成了无源之水。
2. **Bessel 方程里的 $\tau$ 与圆盘问题里的 $r$ 不是同一个变量，讲稿混用了。** 讲稿先写 $v$ 满足的方程（自变量 $r$、含 $\lambda$），紧接着写「$f''+\frac1\tau f'+(1-\frac{n^2}{\tau^2})f=0$」并称之为 Bessel 方程。这两式之间差一个**伸缩换元** $\tau=\sqrt\lambda\,r$（不是把 $r$ 改名成 $\tau$，也不是恒等变形）。不点明这一点，会以为 $\lambda$ 从方程里消失了。
3. **$(1-\frac{n^2}{\tau^2})$ 里的 $n$ 与圆盘问题里的角量子数 $m$ 是同一个量，讲稿用了两个字母。** 本篇统一写 $n=m\in\mathbb Z_{\ge0}$，以免误以为出现了新的整数指标。
4. **「这是分离变量法的理论基础」说得不够。** 乘积是特征函数 ≠ 全部特征函数都是乘积，完备性要另证（§17.1 的要点给出了论证）。若不补这一步，「矩形上所有特征函数」这句话就没有依据，「用分离变量解矩形上的热方程」也就只是猜出了一个解。
5. **圆盘分离变量里「$w$ 为三角函数」跳过了 $k=m^2$ 的论证。** 单值性（$2\pi$-周期）是圆盘问题特有的约束，它直接决定 $k$ 只能取平方数；不写这一步，读者无法理解为什么圆盘的特征值要比「任意 $k$」稀疏得多。
6. **矩形部分只写了 Dirichlet 的例子。** 若换边值条件，$\mu$、$\rho$ 要取对应的谱（Neumann 时含 $0$ 与常数 $1/\sqrt a$），此时矩形上会出现 $\lambda=0$ 的常数特征函数与「只在一个方向振动」的模态。这在一维已由 [[15-Sturm-Liouville定理与一维特征值问题]] 处理，二维只需逐方向套用，本篇未重复。

---

系列：[[16-特征值的变分刻画与Poincaré不等式]] · [[15-Sturm-Liouville定理与一维特征值问题]] · [[18-热方程·建模、边值条件与分离变量法]]
