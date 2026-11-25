# 14 谱理论·鼓面的数学建模与 Laplace 特征值

> 系列：[[13-奇点可去性与局部解析性]] · [[15-Sturm-Liouville定理与一维特征值问题]] · [[02-边值问题的开始·三类边界条件与Dirichlet原理]]

到这里为止，本课一直在解**一个**方程（$\Delta u=f$ 或 $\Delta u=0$）。从本篇开始换问题：**不给定右端项，而是问「哪些 $\lambda$ 使 $-\Delta u=\lambda u$ 有非零解」**。这一支叫**谱理论**，它的物理出处就是鼓面。

讲稿的路线是：先由一个具体的力学模型逼出特征值问题，再抽象成定义与谱定理，最后回到模型解释「音调」与「振型」。本篇按这条线走。

---

## 14.1 鼓面的数学建模

设 $\Omega\subseteq\mathbb R^2$ 是一块区域，$\partial\Omega$ 固定（鼓皮绷在鼓框上）。令 $u(x,t)$ 表示 $t$ 时刻 $x$ 处鼓面的**高度**（偏离平衡位置，向上为正）。

**假设（1）（同步振动）**　每一点的振动节奏相同，只是振幅不同：

$$
u(x,t)=\phi(x)h(t).
$$

**假设（2）（均匀膜）**　面密度 $m>0$，各处拉力系数 $\tau>0$（$\tau$ 由材料决定）。

受力分析：$t$ 时刻 $x$ 处的回复力，来自四周拉力的倾斜分量，方向与位移相反、大小正比于位移，记为

$$
F(x,t)=-C\,\phi(x)h(t),\qquad C>0 ,
$$

$C$ 称为**恢复力系数**（它不是一个独立材料常数，最后会由频率定出来，见下）。

### 14.1.1 时间部分

把 $F$ 与牛顿第二定律放在一起（$F_3$ 表示第三方向即法向的分力）：

$$
F_3(x_0,t)=m(x_0)\left(\frac{d}{dt}\right)^{2}u(x_0,t)=m(x_0)\phi(x_0)h''(t),
$$

与 $F=-C\phi h$ 比较得

$$
-C\,\phi(x_0)h(t)=m\,\phi(x_0)h''(t)\quad\Longrightarrow\quad h''(t)=-\frac Cm\,h(t).
$$

以 $\omega:=\sqrt{C/m}$ 为角频率，方程的通解是

$$
h(t)=h(0)\cos\omega t+\frac{h'(0)}{\omega}\sin\omega t .
$$

### 14.1.2 空间部分

任取 $x_0\in\Omega$，考虑以 $x_0$ 为心、边长 $\delta>0$ 的小正方形，取 $e_1=(1,0)$、$e_2=(0,1)$。以 $e_1$ 方向为例，右侧边界上的拉力为

$$
F_{\text{右},\delta}(x_0,t)=\tau\bigl[u(x_0+\delta e_1,t)-u(x_0,t)\bigr],
$$

（拉力正比于相邻两点的**高度差**，这正比于斜率）。四个方向相加，这个小正方形在 $x_0$ 处所受的力为

$$
F_\delta(x_0,t)=\tau h(t)\Bigl[\phi(x_0+\delta e_1)-\phi(x_0)+\phi(x_0-\delta e_1)-\phi(x_0)+\phi(x_0+\delta e_2)-\phi(x_0)+\phi(x_0-\delta e_2)-\phi(x_0)\Bigr].
$$

单位面积受力即 $\lim_{\delta\to0}F_\delta/\delta^2$，而方括号除以 $\delta^2$ 正是 $\Delta\phi(x_0)$ 的二阶中心差分：

$$
\lim_{\delta\to0}\frac{F_\delta(x_0,t)}{\delta^2}=\tau h(t)\,\Delta\phi(x_0).
$$

### 14.1.3 合成

时间部分给 $F=-C\phi h$，空间部分给 $F=\tau h\Delta\phi$。两者相等，约去 $h(t)\phi(x_0)$：

$$
-C\,\phi(x_0)h(t)=\tau h(t)\Delta\phi(x_0)\quad\Longrightarrow\quad-\Delta\phi=\frac{C}{\tau}\,\phi .
$$

合起来（$\lambda:=C/\tau$）：

$$
\begin{cases}h''(t)+\dfrac Cm h(t)=0,\\[6pt] -\Delta\phi=\lambda\phi, & x\in\Omega,\\[4pt] \phi|_{\partial\Omega}\equiv0,\end{cases}\tag{14.1}
$$

最后一条是「鼓框固定」：边界上高度恒为零，即 **Dirichlet 边值条件**。

角频率与频率：

$$
\omega=\sqrt{\frac Cm}=\sqrt{\frac{\tau\lambda}{m}},\qquad f=\frac{\omega}{2\pi}.
$$

**频率决定音调**（$f$ 大则音高），**振型 $\phi$ 决定这个音是怎么振的**（节线在哪）。

### 14.1.4 「是不是每个 $\phi$ 都能配一个 $\lambda$」——不能

一般情形下，真正的方程是

$$
m\left(\frac{d}{dt}\right)^{2}u(x,t)=\tau\,\Delta_x u(x,t),
$$

它对**任意**初值都成立；但给定 $u(x,0)$ 之后，$u(x,t)$ 一般**不能**分离成 $\phi(x)h(t)$。能分离的只有形如 $\phi(x)h(t)$ 的**简单振动**（驻波），而一般的振动是它们的叠加：

$$
u(x,t)=\sum_{i=1}^{\infty}\phi_i(x)h_i(t).
$$

也就是说：**「每个 $\phi$ 都有一个 $\lambda$」是错的，正确的是「只有一列特殊的 $\phi$ 才配得上 $\lambda$」**。反过来，$-\Delta$ 在 Dirichlet 条件下恰好给出一列（可数）特征函数，它们构成 $L^2(\Omega)$ 的基 —— 这就是后面谱定理的内容，也正是「叠加为什么合法」的依据。

**听音辨鼓**：知道鼓面产生**全部**音调（即知道整个谱 $\lambda(\Omega)$），能不能还原鼓的形状 $\Omega$？这是 Mark Kac 1966 年提出的著名问题（"Can one hear the shape of a drum?"），答案是否定的（Gordon–Webb–Wolpert, 1992 给出了等谱非等距的区域）。

---

## 14.2 特征值、特征函数、特征对

**定义 14.2.1**　称 $\lambda$ 为区域 $\Omega\subseteq\mathbb R^n$ 上 **Dirichlet 边值条件**下的特征值，若存在非零函数 $\phi\in C^2(\Omega)\cap C(\overline\Omega)$ 满足

$$
\begin{cases}-\Delta\phi=\lambda\phi, & x\in\Omega,\\ \phi=0, & x\in\partial\Omega.\end{cases}
$$

**Neumann 边值条件**：

$$
\begin{cases}-\Delta\phi=\lambda\phi, & x\in\Omega,\\[4pt] \dfrac{\partial\phi}{\partial\nu}=0, & x\in\partial\Omega.\end{cases}
$$

**Robin 边值条件**（$\sigma>0$）：

$$
\begin{cases}-\Delta\phi=\lambda\phi, & x\in\Omega,\\[4pt] \dfrac{\partial\phi}{\partial\nu}+\sigma\phi=0, & x\in\partial\Omega.\end{cases}
$$

称相应的 $\phi$ 为**特征函数**，$(\lambda,\phi)$ 为 $\Omega$ 上的**特征对**。

> **要点（三件事必须记住）**
> 1. **$\phi\not\equiv0$ 是定义的一部分。** 少了这一条，「任意 $\lambda$ 都是特征值」（因为 $\phi\equiv0$ 永远满足方程）。
> 2. **涉及 $\partial\phi/\partial\nu$ 的边值条件要求特征函数有相应的边界正则性。** 上面的 $C^2(\Omega)\cap C(\overline\Omega)$ 对 Neumann、Robin 是不够的（要 $\nabla\phi$ 连续到边界），严格写法要用 $H^2(\Omega)$；课堂在有界光滑区域上陈述，本篇沿用。
> 3. **三类条件的差别全在边界项上**，这一点在 §14.3 的能量恒等式里看得最清楚。

**注（特征空间）**

1. 若 $\phi$ 是 $\lambda$ 的特征函数，则非零倍数 $a\phi$ 也是；若 $\phi_1,\phi_2$ 都对应 $\lambda$，则 $a\phi_1+b\phi_2$ 仍满足同一方程与边值条件。把零函数一并加进去，所得线性空间 $E_\lambda:=\{u\in L^2(\Omega):-\Delta u=\lambda u,\ u\ \text{满足所选边值条件}\}$ 称为**特征空间**。
2. **类比矩阵**：对称 $n$ 阶矩阵 $A$ 的特征值只有 $n$ 个（记 $\lambda(A)=\{\lambda_1,\dots,\lambda_n\}$），且特征向量可正交对角化。下面的谱定理说的正是「$-\Delta$ 是一个有无穷多个特征值、且有完备正交特征函数的『对称算子』」——算子的对称性由 Green 第一恒等式（[[08-Green恒等式与Green表示公式]]）保证，而**对称**正是特征值全为实、特征函数相互正交的根源。

---

## 14.3 特征值的符号

**定理 14.3.1**　Dirichlet 特征值满足 $\lambda>0$；Neumann 特征值满足 $\lambda\ge0$。若 $\Omega$ 连通，则 Neumann 的零特征函数为非零常数。

**证明**　由 $-\Delta\phi=\lambda\phi$，两边同乘 $\phi$ 并在 $\Omega$ 上积分：

$$
\lambda\int_\Omega\phi^2\,dx=\int_\Omega-\phi\Delta\phi\,dx .
$$

把被积函数拆成「散度项 + 平方项」。由 $\operatorname{div}(\phi\nabla\phi)=\phi\Delta\phi+|\nabla\phi|^2$，即 $-\phi\Delta\phi=-\operatorname{div}(\phi\nabla\phi)+|\nabla\phi|^2$，逐项写开：

$$
\lambda\int_\Omega\phi^2\,dx=-\sum_{j=1}^{n}\int_\Omega(\phi\phi_{x_j})_{x_j}\,dx+\int_\Omega|\nabla\phi|^2\,dx=-\int_{\partial\Omega}\phi\,\frac{\partial\phi}{\partial\nu}\,d\sigma+\int_\Omega|\nabla\phi|^2\,dx .
$$

（第二步用散度定理，$\sum_j\phi\phi_{x_j}\nu_j=\phi\,\partial_\nu\phi$。）于是

$$
\lambda\int_\Omega\phi^2\,dx=\int_\Omega|\nabla\phi|^2\,dx-\int_{\partial\Omega}\phi\,\frac{\partial\phi}{\partial\nu}\,d\sigma .
$$

**Dirichlet**：$\phi|_{\partial\Omega}=0$，边界项为零，故 $\lambda\int\phi^2=\int|\nabla\phi|^2\ge0$。若 $\lambda=0$，则 $\nabla\phi\equiv0$，$\phi$ 为常数，而 $\phi|_{\partial\Omega}=0$ 迫使 $\phi\equiv0$，与 $\phi\not\equiv0$ 矛盾。故 $\lambda>0$。

**Neumann**：$\partial_\nu\phi|_{\partial\Omega}=0$，边界项同样为零，$\lambda\ge0$。若 $\lambda=0$，则 $\nabla\phi\equiv0$，即 $\phi$ 为常数，而常数确实满足 $\partial_\nu\phi=0$ —— 只要 $\Omega$ 连通，常数就是全部解，$\phi\equiv$ 某非零常数。$\square$

> **要点**　这条证明是**整个谱理论唯一需要的技术**：把 $\lambda\int\phi^2$ 转化成 $\int|\nabla\phi|^2$ 与一个边界项。后面 [[16-特征值的变分刻画与Poincaré不等式]] 的 Rayleigh 商、Robin 条件下的 $\lambda>0$、以及 Poincaré 不等式的最佳常数，全是这一式的变形。
>
> **补充（Robin 情形）**　定理只处理了 Dirichlet 与 Neumann。Robin 条件的能量计算见 [[16-特征值的变分刻画与Poincaré不等式]] §16.1，结论是 $\lambda>0$（$\sigma>0$ 时边界项 $-\int\phi\partial_\nu\phi=\sigma\int_{\partial\Omega}\phi^2>0$，故不仅非负而且严格正）。

---

## 14.4 $L^2(\Omega)$ 与完备标准正交基

**定义 14.4.1**　

$$
L^2(\Omega):=\left\{f:\Omega\to\mathbb R\ \text{可测}\ \Big|\ \int_\Omega|f|^2\,dx<+\infty\right\},
$$

其中几乎处处相等的函数视为同一元素。

**定义 14.4.2**　称 $\{f_j\}_{j=1}^\infty$ 为 $L^2(\Omega)$ 的一组**完备标准正交基**，若 $\{f_j\}\subseteq L^2(\Omega)$ 满足

1. $\displaystyle\int_\Omega f_if_j\,dx=\delta_{ij}$（标准且正交）；
2. $\overline{\operatorname{Span}\{f_j\}_{j=1}^\infty}=L^2(\Omega)$，闭包取在 $L^2(\Omega)$ 拓扑下。

> **要点**　条件 2 是「完备」二字的全部内容：**不是**说每个 $L^2$ 函数都是有限线性组合，而是说有限线性组合在 $L^2$ 范数下稠密。所以展开式 $\sum a_jf_j=f$ 只在 $L^2$ 意义下成立，**不一定逐点收敛**，更不一定连续。这一点在 [[18-热方程·建模、边值条件与分离变量法]] 讲一维分离变量时会被反复用到。

---

## 14.5 谱定理

**定理 14.5.1（Laplace 算子的谱）**　设 $\Omega\subseteq\mathbb R^n$ 为有界光滑区域。带上述三类边值条件之一的 Laplace 算子，其特征值（按重数排列）构成离散集

$$
0\le\lambda_1\le\lambda_2\le\cdots\le\lambda_k\le\cdots\longrightarrow+\infty ,
$$

且相应的特征函数满足：

1. $\{u_k\}_{k=1}^\infty$ 可构成 $L^2(\Omega)$ 的一组完备标准正交基；
2. 每个特征空间维数有限：

$$
\dim E_\lambda<\infty,\qquad E_\lambda:=\{u\in L^2(\Omega):-\Delta u=\lambda u,\ u\ \text{满足所选边值条件}\}.
$$

其中 $\lambda_1=0$ 只在 Neumann 条件（$\Omega$ 连通）下出现，此时特征函数为非零常数；Dirichlet 与 Robin（$\sigma>0$）都有 $\lambda_1>0$。

**补充（为什么会有离散谱 —— 讲稿只说「用紧算子性质，自行查阅」）**

把 $-\Delta$ 的逆算子写出来：给定 $f\in L^2(\Omega)$，解

$$
\begin{cases}-\Delta u=f, & x\in\Omega,\\ u=0, & x\in\partial\Omega,\end{cases}
$$

并记 $Tf:=u$。则：

1. **$T$ 是紧算子**：椭圆正则性给出 $\lVert u\rVert_{H^2}\le C\lVert f\rVert_{L^2}$，而 $H^2(\Omega)\hookrightarrow L^2(\Omega)$ 是紧嵌入（Rellich–Kondrachov），故 $T:L^2\to L^2$ 紧；
2. **$T$ 是自伴的**：由 Green 第一恒等式，$\int(Tf)g=\int f(Tg)$；
3. **$T$ 是正算子**：$\int(Tf)f=\int|\nabla Tf|^2\ge0$。

紧自伴算子的谱定理给出：$T$ 的特征值 $\mu_k\searrow0$（$\mu_k>0$）、每个特征空间有限维、特征向量构成 $L^2$ 的完备正交基。把 $\mu_k=1/\lambda_k$ 翻译回来，就是上面定理的三条结论。**注意这一路只对 Dirichlet 直接成立**；Neumann / Robin 要用 $-\Delta+1$ 或商空间（常数函数之外）来避开核。

这条补充只交代机制，细节（$H^2$ 正则性、嵌入定理）超出本课范围。

**在鼓面模型中**：称

$$
\lambda(\Omega):=\{0<\lambda_1\le\lambda_2\le\cdots\le\lambda_k\le\cdots\to+\infty\}
$$

为 $\Omega$ 上 Laplace 算子的**谱**（鼓面是 Dirichlet 条件，故从正数开始）。

**听音辨鼓问题**：已知 $\lambda(\Omega)$，能否还原 $\Omega$ 的形状？（答案是否定的，见 §14.1.4。）

---

## 14.6 回到鼓面：模态展开

对任意敲击方式产生 $u(x,t)$。一方面牛顿第二定律给

$$
\left(\frac{d}{dt}\right)^2u(x,t)=\frac1mF(x,t);
$$

另一方面受力分析给

$$
\Delta_xu(x,t)=\frac1\tau F(x,t).
$$

合并（消去 $F$）得

$$
\left(\frac{\partial}{\partial t}\right)^2u(x,t)=\Delta_xu(x,t)\cdot\frac{\tau}{m}.
$$

取谱定理给的标准正交特征函数 $\phi_k$。对固定时刻 $t$，把 $u(\cdot,t)$ 按 $\{\phi_k\}$ 展开：

$$
u(x,t)=\sum_{k=1}^\infty q_k(t)\phi_k(x),
$$

级数在 $L^2(\Omega)$ 中收敛，系数由正交性给出

$$
q_k(t)=\int_\Omega u(x,t)\phi_k(x)\,dx .
$$

对可以求导两次的（光滑）解，代入方程并比较系数：

$$
\begin{cases}q_k''(t)=\dfrac{\tau}{m}\,(-\lambda_k)\,q_k(t),\\[6pt] -\Delta\phi_k=\lambda_k\phi_k .\end{cases}
$$

也就是每个模态独立地做简谐振动，角频率 $\omega_k=\sqrt{\tau\lambda_k/m}$ —— 与 §14.1.3 的 $\omega=\sqrt{\tau\lambda/m}$ 完全一致。

> **要点**　这就是「鼓的音调 = $\sqrt{\tau\lambda_k/m}/(2\pi)$」的严格版本：**谱 $\{\lambda_k\}$ 决定全部音高，特征函数 $\{\phi_k\}$ 决定每个音的振型（节线位置）**。而 $\tau/m$ 只是一个整体缩放因子，改变它只会把所有音高同比例升降（就像把鼓皮绷紧），不改变音程关系 —— 所以「听音辨鼓」问的确实是 $\lambda(\Omega)$ 与 $\Omega$ 的关系。

---

## 订正

1. **时间方程里把恢复力系数 $C$ 写成了张力 $\tau$。** 讲稿 §14.1.1 与 §14.1.3 共有四处把第一个 $\tau$ 当成回复力系数：写作 $-\tau\phi(x_0)h(t)=m\phi(x_0)h''(t)$、$h''(t)=-\dfrac{\tau}{m}h(t)$、$h(t)=h(0)\cos\bigl(\sqrt{\tau/m}\,t\bigr)+h'(0)\sqrt{m/\tau}\sin\bigl(\sqrt{\tau/m}\,t\bigr)$，以及 $-\tau\phi h=\tau h\Delta\phi$。四处都应把第一个 $\tau$ 换成 $C$，即 $F=-C\phi h$、$h''=-\dfrac Cm h$、$h(t)=h(0)\cos\omega t+\dfrac{h'(0)}{\omega}\sin\omega t$（$\omega=\sqrt{C/m}$）、$-C\phi h=\tau h\Delta\phi$。
**判据**：讲稿自己在 (14.1) 里写的是 $h''+\frac Cmh=0$ 与 $-\Delta\phi=\frac C\tau\phi$，与上面的 $\tau$ 不自洽；若真取 $C=\tau$，则 $\lambda\equiv1$，所有模态频率相同，与「鼓有好几个音」直接矛盾。**$C$ 不是材料常数，它等于 $m\omega^2$，是由频率反定的量。**
2. **谱定理括注里的「Robin」应为「Neumann」。** 讲稿写「Neumann 条件容许 $\lambda_1=0$；这里的 Robin 条件取 $\sigma>0$」，前半句是对的、后半句放在「容许 $\lambda_1=0$」这个位置上就错位了 —— Robin 且 $\sigma>0$ 时 $\lambda_1>0$（见 [[16-特征值的变分刻画与Poincaré不等式]] §16.1）。这句话本身在讲稿里也带了自我更正（「手稿括注中的 Robin 应为 Neumann」），本篇按正确版写。
3. **定理 14.5.1 的 $\lambda_1=0$ 判别少了「$\Omega$ 连通」。** 不连通时 Neumann 的零特征空间维数等于 $\Omega$ 的连通分支数，$\lambda_1=0$ 仍成立，但零特征函数是「分片常数」而非常数。定理正文已补。
4. **定义 14.2.1 只写 $\phi\in C^2(\Omega)\cap C(\overline\Omega)$ 对 Neumann / Robin 不够。** 边界条件里出现 $\partial_\nu\phi$，要求 $\nabla\phi$ 连续到 $\partial\Omega$；严格陈述应在 $H^2(\Omega)$ 中做。这是讲稿自己也点出的问题（「涉及法向导数时要求相应边界正则性」），此处仅记录。
5. **『手稿将 $\phi$ 直观地称为「初始时刻敲击鼓面的力度」』这个说法是错的。** $\phi$ 是**振型**（空间模态），由 $\Omega$ 与 $\lambda$ 唯一决定，与怎么敲无关；「敲击的力度」体现在展开系数 $q_k(0)$ 与 $q_k'(0)$ 上（即初始位移与初始速度）。同一个 $\phi$ 配不同的 $q_k(0)$，得到的是同一个音、不同响度 —— 响度对 $\phi$ 没有影响。
6. **§14.1.2 里 $F_\delta/\delta^2\to\tau h\Delta\phi$ 的「单位面积」没说清。** 严格地说，$F_\delta$ 是小正方形四条边上拉力的合力（量纲为力），除以面积 $\delta^2$ 才是面力密度；而 $\Delta\phi$ 是 $\delta\to0$ 时二阶中心差分的极限，两者差一个 $\tau h$ 的因子，正是拉力系数的量纲来源。

---

系列：[[13-奇点可去性与局部解析性]] · [[15-Sturm-Liouville定理与一维特征值问题]] · [[02-边值问题的开始·三类边界条件与Dirichlet原理]]
