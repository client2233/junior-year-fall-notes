# 08 Green 恒等式与 Green 表示公式

> 系列：[[07-基本解与广义函数]] · [[09-Green函数的性质]] · [[10-球与半空间上的Green函数]]

本篇回答一个问题：**给定 $\Omega$ 上有界、$\partial\Omega$ 为 $C^1$ 曲面，Dirichlet 问题是否存在解？**

路线是「先造一个公式」。从 Green 恒等式出发，把区域内的 $u(x)$ 用源项 $-\Delta u$ 与**边界上**的 $u$、$\partial u/\partial\nu$ 表示出来——这就是 Green 表示公式。它本身还不解决存在性（因为它同时用到了两类边界数据，条件过剩），但它指出了下一步该构造什么。

---

## 8.1 Green 恒等式

**定理（Green 恒等式）**　令 $\Omega\subseteq\mathbb R^n$ 为有界、边界 $C^1$ 的区域，记 $\nu$ 为单位外法向量，$u,v\in C^2(\Omega)\cap C^1(\overline\Omega)$。则：

**（1）Green 第一公式**

$$
\int_\Omega u\cdot\Delta v\,dx=-\int_\Omega\langle\nabla u,\nabla v\rangle\,dx+\int_{\partial\Omega}u\frac{\partial v}{\partial\nu}\,d\sigma_x .
$$

**证明**

$$
\int_\Omega u\cdot\Delta v\,dx
=\sum_{i=1}^n\int_\Omega u\cdot v_{x_ix_i}\,dx
=\sum_{i=1}^n\int_\Omega (u\cdot v_{x_i})_{x_i}\,dx-\int_\Omega\langle\nabla u,\nabla v\rangle\,dx
=\int_{\partial\Omega}u\cdot\langle\nabla v,\nu\rangle\,d\sigma_x-\int_\Omega\langle\nabla u,\nabla v\rangle\,dx
=\int_{\partial\Omega}u\frac{\partial v}{\partial\nu}\,d\sigma_x-\int_\Omega\langle\nabla u,\nabla v\rangle\,dx .
$$

（第三、四个等号分别用了散度定理 $\int_\Omega\operatorname{div}F=\int_{\partial\Omega}\langle F,\nu\rangle$ 与 $\dfrac{\partial v}{\partial\nu}=\langle\nabla v,\nu\rangle$。）$\square$

**（2）Green 第二公式**

$$
\int_\Omega\bigl(u\Delta v-v\Delta u\bigr)dx
=\int_{\partial\Omega}\left(u\frac{\partial v}{\partial\nu}-v\frac{\partial u}{\partial\nu}\right)d\sigma_x .
$$

> **补充（第二公式的来历）**　把第一公式中 $u,v$ 互换，得到

$$
\int_\Omega v\cdot\Delta u\,dx=-\int_\Omega\langle\nabla v,\nabla u\rangle\,dx+\int_{\partial\Omega}v\frac{\partial u}{\partial\nu}\,d\sigma_x ;
$$

> 用原式减去它，两侧的 $\int\langle\nabla u,\nabla v\rangle$ 完全相同而相消，即得第二公式。第二公式的特点是**左端关于 $u,v$ 反对称、右端只含边界**——这正是它能被用来「反解出 $u(x)$」的原因。

> **注（允许 $\Delta u$ 无界）**　讲稿明确指出：**允许 $\Delta u$ 在趋近边界时无界（包括趋于 $\infty$），也不要求 $\Delta u\in L^1(\Omega)$**，对 $\Delta v$ 同理。这里及下文含 $\Delta u$ 或 $\Delta v$ 的体积分均按**广义积分**理解：取 $\Omega_j\Subset\Omega$、$\Omega_j\uparrow\Omega$，其光滑边界在 $C^1$ 意义下逼近 $\partial\Omega$，相应积分定义为在 $\Omega_j$ 上积分的极限；基本解的奇点按**去球极限**处理。

---

## 8.2 Green 表示公式

**回忆基本解**（见 [[07-基本解与广义函数]]）：

$$
\Gamma(x):=\begin{cases}
-\dfrac{1}{2\pi}\log|x|, & n=2,\\[8pt]
\dfrac{1}{(n-2)|S^{n-1}|\,|x|^{n-2}}, & n\ge3,
\end{cases}
\qquad
-\Delta_x\Gamma(x)=\delta_0,\qquad -\Delta_y\Gamma(x-y)=\delta_x .
$$

**定理（Green 表示公式）**　令 $\Omega\subseteq\mathbb R^n$ 为有界、边界 $C^1$ 的区域，$u\in C^2(\Omega)\cap C^1(\overline\Omega)$，则对于 $\forall x\in\Omega$，

$$
u(x)=\int_\Omega\bigl[-\Delta u(y)\bigr]\Gamma(x-y)\,dy
+\int_{\partial\Omega}\left[\Gamma(x-y)\frac{\partial u}{\partial\nu_y}(y)-u(y)\frac{\partial\Gamma(x-y)}{\partial\nu_y}\right]d\sigma_y .
$$

**证明**　对于 $\forall x\in\Omega$，取 $B_r(x)\Subset\Omega$，记 $\Omega_r:=\Omega\setminus\overline{B_r(x)}$。在 $\Omega_r$ 上对 $u(y)$ 与 $v(y)=\Gamma(x-y)$ 用 Green 第二公式。由于 $y\ne x$ 于 $\Omega_r$，有 $\Delta_y\Gamma(x-y)=0$，故

$$
\int_{\Omega_r}\bigl[-\Delta u(y)\bigr]\Gamma(x-y)\,dy
=\int_{\Omega_r}\bigl[-\Delta u(y)\bigr]\Gamma(x-y)\,dy+\int_{\Omega_r}u(y)\cdot\Delta_y\Gamma(x-y)\,dy .
$$

（第二项为零，写出来是为了凑成 $u\Delta v-v\Delta u$ 的形式。）于是

$$
=\int_{\partial\Omega}\left[u(y)\frac{\partial\Gamma(x-y)}{\partial\nu_y}-\Gamma(x-y)\frac{\partial u}{\partial\nu_y}(y)\right]d\sigma_y
-\int_{\partial B_r(x)}\left[u(y)\frac{\partial\Gamma(x-y)}{\partial\nu_y}-\Gamma(x-y)\frac{\partial u}{\partial\nu_y}(y)\right]d\sigma_y .
$$

> 第二项前的负号来自：在 $\partial B_r(x)$ 上，作为 $\Omega_r$ 的外法向量有 $\nu_r=-\nu_B$（$\nu_B$ 为小球的单位外法向量，即由 $x$ 指向球外），故 $\partial\Gamma/\partial\nu_r=-\partial\Gamma/\partial\nu_B$。

**把奇点附近的项算出来。** 记 $r=|x-y|$，在 $\partial B_r(x)$ 上 $\Gamma(x-y)=\Gamma(r)$ 是常数，且（用 $\nabla_y\Gamma(x-y)=-\nu_B\big/\bigl(|S^{n-1}|r^{n-1}\bigr)$）

$$
\frac{\partial\Gamma}{\partial\nu_y}(x-y)=\begin{cases}
-\dfrac{1}{2\pi r}, & n=2,\\[8pt]
-\dfrac{1}{|S^{n-1}|r^{n-1}}, & n\ge3,
\end{cases}
\qquad\text{即}\qquad
\frac{\partial\Gamma}{\partial\nu_y}\Big|_{\partial B_r(x)}=-\frac{1}{|\partial B_r(x)|}.
$$

代入得

$$
\int_{\Omega_r}\bigl[-\Delta u(y)\bigr]\Gamma(x-y)\,dy=:\ I+II+III,
$$

其中

$$
I=\int_{\partial\Omega}\left[u(y)\frac{\partial\Gamma(x-y)}{\partial\nu_y}-\Gamma(x-y)\frac{\partial u}{\partial\nu_y}(y)\right]d\sigma_y,
II=\frac{1}{|\partial B_r(x)|}\int_{\partial B_r(x)}u(y)\,d\sigma_y,\qquad
III=\Gamma(r)\int_{\partial B_r(x)}\frac{\partial u}{\partial\nu_y}\,d\sigma_y .
$$

**取极限。** 由 $u$ 在 $x$ 处的连续性，$II\to u(x)$；又由 $u\in C^1(\overline\Omega)$ 与 $\Gamma(r)|\partial B_r(x)|\to0$，

$$
|III|\le\Gamma(r)\cdot|\partial B_r(x)|\cdot\sup_{\partial B_r(x)}|\nabla u|\longrightarrow0 .
$$

于是令 $r\to0^+$，此时左端趋于 $\displaystyle\int_\Omega[-\Delta u(y)]\Gamma(x-y)dy$，得

$$
\int_\Omega\bigl[-\Delta u(y)\bigr]\Gamma(x-y)\,dy
=\int_{\partial\Omega}\left[u\frac{\partial\Gamma}{\partial\nu_y}-\Gamma\frac{\partial u}{\partial\nu_y}\right]d\sigma_y+u(x),
$$

即

$$
u(x)=\int_\Omega\bigl[-\Delta u(y)\bigr]\Gamma(x-y)\,dy
+\int_{\partial\Omega}\left[\Gamma(x-y)\frac{\partial u}{\partial\nu_y}(y)-u(y)\frac{\partial\Gamma(x-y)}{\partial\nu_y}\right]d\sigma_y .
$$

$\square$

> **要点**　表示公式的读法是：**区域内一点的函数值 = 源的贡献（体积分）+ 边界上的「单层势」$\Gamma\partial u/\partial\nu$ 与「双层势」$u\partial\Gamma/\partial\nu$**。它的推导只做了一件事——把 $\Delta$ 从 $u$ 挪到已知核 $\Gamma$ 上，再用挖球把 $\Gamma$ 的奇性单独处理掉。这与 [[01-从万有引力到Poisson方程]] 定理 1.1 的证明是同一个套路。

> **补充（调和函数的复制公式）**　若在 $\Omega$ 内 $\Delta u=0$，表示公式退化为

$$
u(x)=\int_{\partial\Omega}\left[\Gamma(x-y)\frac{\partial u}{\partial\nu_y}(y)-u(y)\frac{\partial\Gamma(x-y)}{\partial\nu_y}\right]d\sigma_y\qquad(x\in\Omega),
$$

> 即**调和函数完全由它在边界上的值与法向导数决定**。这条公式后面在球上就变成 Poisson 积分公式（[[10-球与半空间上的Green函数]]、[[11-Poisson积分公式与Dirichlet问题的可解性]]）。

---

## 8.3 边界条件的「过剩」与 Green 函数的引出

**观察**　令 $\Omega\subseteq\mathbb R^n$，$\partial\Omega$ 为 $C^1$ 曲面且 $\Omega$ 有界。若存在 $u\in C^2(\Omega)\cap C^1(\overline\Omega)$ 满足 Dirichlet 问题

$$
\begin{cases}
-\Delta u=f, & x\in\Omega,\\
u=g, & x\in\partial\Omega,
\end{cases}
\qquad\text{或}\qquad
\begin{cases}
-\Delta u=f, & x\in\Omega,\\[4pt]
\dfrac{\partial u}{\partial\nu}=h, & x\in\partial\Omega,
\end{cases}
$$

则 $u$ 满足（把表示公式中的 $-\Delta u$ 写成 $f$、$u|_{\partial\Omega}$ 写成 $g$、$\partial u/\partial\nu|_{\partial\Omega}$ 写成 $h$）

$$
u(x)=\int_\Omega f(y)\Gamma(x-y)\,dy
+\int_{\partial\Omega}\left[-g(y)\frac{\partial\Gamma(x-y)}{\partial\nu_y}+h(y)\Gamma(x-y)\right]d\sigma_y .
$$

**但这个公式对 Dirichlet 问题的存在性是「过剩」的**：解的函数值由 $-\Delta u$、$u|_{\partial\Omega}$、$\left.\dfrac{\partial u}{\partial\nu}\right|_{\partial\Omega}$ **三个**函数的积分共同决定。而 [[02-边值问题的开始·三类边界条件与Dirichlet原理]] 与 [[05-下调和函数·弱极值原理与Hopf边界点引理]] 告诉我们，Dirichlet 问题与 Neumann 问题各自都有唯一性（后者在相差常数的意义下）。**既然解由两类边值条件中的任一类就能唯一确定，那么把它写成同时依赖两类条件的积分，条件就是多余的。** 另外，两类边值条件也不一定相容（例如 Neumann 问题需要 $\int_\Omega f=\int_{\partial\Omega}g$，见 [[05-下调和函数·弱极值原理与Hopf边界点引理]] 的补充）。

### 解 Dirichlet 问题：消掉 $\partial u/\partial\nu$

既然 Dirichlet 问题的边界数据是 $u|_{\partial\Omega}=g$，我们希望表示公式里**不要出现** $\dfrac{\partial u}{\partial\nu}$。办法是选一个「修正项」把 $\Gamma$ 改掉：

> **若能找到函数 $\gamma(x,y)$ 满足：对 $\forall x\in\Omega$，**

$$
\begin{cases}
\Delta_y\gamma(x,y)=0, & y\in\Omega,\\[4pt]
\gamma(x,y)=\Gamma(x-y), & y\in\partial\Omega,
\end{cases}
$$

> **则令**

$$
G(x,y):=\Gamma(x-y)-\gamma(x,y),
$$

> **就可以把 $\partial u/\partial\nu$ 从公式中消掉。**

**推导**　用表示公式写出 $u(x)$，再对 $u$ 与 $\gamma$ 用 Green 第二公式（此时两者都在 $\Omega$ 上光滑，无需挖球）：

$$
\int_\Omega\bigl[-\Delta u(y)\bigr]G(x,y)\,dy
=\int_\Omega\bigl[-\Delta u(y)\bigr]\Gamma(x-y)\,dy+\int_\Omega\Delta u(y)\gamma(x,y)\,dy
=u(x)+\int_{\partial\Omega}\left[u(y)\frac{\partial\Gamma(x-y)}{\partial\nu_y}-\Gamma(x-y)\frac{\partial u}{\partial\nu_y}(y)\right]d\sigma_y
+\int_{\partial\Omega}\left[\gamma(x,y)\frac{\partial u}{\partial\nu_y}(y)-u(y)\frac{\partial\gamma(x,y)}{\partial\nu_y}\right]d\sigma_y
=u(x)+\int_{\partial\Omega}\left[u(y)\frac{\partial G(x,y)}{\partial\nu_y}-G(x,y)\frac{\partial u}{\partial\nu_y}(y)\right]d\sigma_y .
$$

而在 $\partial\Omega$ 上 $G(x,y)=\Gamma(x-y)-\gamma(x,y)=0$，故含 $\dfrac{\partial u}{\partial\nu}$ 的项整体消失：

$$
\int_\Omega\bigl[-\Delta u(y)\bigr]G(x,y)\,dy=u(x)+\int_{\partial\Omega}u(y)\frac{\partial G(x,y)}{\partial\nu_y}\,d\sigma_y .
$$

即

$$
u(x)=\int_\Omega\bigl[-\Delta u(y)\bigr]G(x,y)\,dy-\int_{\partial\Omega}u(y)\frac{\partial G(x,y)}{\partial\nu_y}\,d\sigma_y .\tag{8.1}
$$

**这就是 Dirichlet 问题的求解公式。**

**分析**　公式 (8.1) 中出现的 $G(x,y)=\Gamma(x-y)-\gamma(x,y)$ **只与 $\Omega$ 有关**——确切说，只与 $\Omega$ 上的 Dirichlet 问题

$$
\begin{cases}
\Delta_y\gamma(x,y)=0, & y\in\Omega,\\[4pt]
\gamma(x,y)=\Gamma(x-y), & y\in\partial\Omega
\end{cases}
$$

**是否可解**有关。只要 $\gamma$ 存在，且 $\gamma\in C^2(\Omega)\cap C^1(\overline\Omega)$，我们便可直接写出 $u(x)$。

**并且——边界正则性**：把 $-\Delta u=f$、$u|_{\partial\Omega}=g$ 代入 (8.1) 并对 $x$ 求法向导数，

$$
\frac{\partial u}{\partial\nu}(x)=\int_\Omega f(y)\frac{\partial G(x,y)}{\partial\nu_x}\,dy-\int_{\partial\Omega}g(y)\frac{\partial^2G(x,y)}{\partial\nu_x\partial\nu_y}\,d\sigma_y .
$$

也即：**$u$ 是否 $C^1(\overline\Omega)$ 取决于 $\gamma$ 是否 $C^2(\overline\Omega)$。**

> **要点**　这一小节的逻辑是整章的关键：**把「解 Dirichlet 问题」归结为「造出 $G$」**，而造 $G$ 归结为「解一个特殊边值数据的 Dirichlet 问题（边界值取 $\Gamma(x-y)$）」。后者的好处是：边界数据 $\Gamma(x-\cdot)$ 是**显式已知**的光滑函数，而原问题的 $g$ 是任意给定的。这就是 [[10-球与半空间上的Green函数]] 只对球和半空间能显式做出 $G$ 的原因——只有在对称性足够高时，$\gamma$ 才能用「对称点」写出来。

> **补充（为什么 $G|_{\partial\Omega}=0$ 是关键）**　注意消项用的是 $G$ 在边界上为零，而不是 $\gamma$ 的某个性质。所以 $G$ 的定义可以理解为：「找一个在 $\partial\Omega$ 上取零的修正核，使得它在 $y\to x$ 时和 $\Gamma(x-y)$ 只差一个有界项」。事实上（见 [[09-Green函数的性质]]）

$$
\lim_{y\to x}\frac{G(x,y)}{\Gamma(x-y)}=1,
$$

> 即 $G$ 保留了 $\Gamma$ 的奇性，同时又在边界上为零——这正是「Green 函数」这个名字的含义。

> **补充（$G$ 的对称性预告）**　(8.1) 里的 $G(x,y)$ 看起来对 $x,y$ 不对称（$x$ 是「源点」、$y$ 是积分变量），但其实有 $G(x,y)=G(y,x)$，见 [[09-Green函数的性质]] §9.3。这条对称性是把 (8.1) 反过来用于 $y$（把 $G(y,x)$ 看作 $x$ 的函数）的关键。

---

## 订正

1. **Green 表示公式的证明漏了最后一步「组装」。** 讲稿的 Pf 以「$\text{原积分}=I+II+III$，其中 $\lim II=u(x)$、$\lim III=0$」结束，**没有写出**由 $\lim(\text{原积分})=I+u(x)$ 反解出

$$
u(x)=\int_\Omega\bigl[-\Delta u\bigr]\Gamma\,dy-I
$$

> 这一步。本篇已在证明中补出。（缺少这一步时，读者容易误以为 $u(x)$ 就是「原积分」本身。）
2. **Green 第二公式讲稿未给推导。** 讲稿直接列出 (1)(2) 两式。第二公式由第一公式互换 $u,v$ 后相减得到，关键点是左侧的 $\int\langle\nabla u,\nabla v\rangle$ 相同而相消。本篇已在**补充**中补出。
3. **定理的陈述应补上对奇性的说明。** 讲稿的定理只写 $u\in C^2(\Omega)\cap C^1(\overline\Omega)$，而公式右端的体积分含 $\Gamma$ 的奇点、边界积分含 $\partial\Gamma/\partial\nu$ 的奇性。讲稿把「去球极限 / 广义积分」的约定写在定理**之前**的注里，本篇保留该约定并说明两处奇性各自由谁处理：体积分的奇点由去球极限处理，边界积分的奇点在 $x\in\Omega$ 固定时其实不奇（$|x-y|>0$）。
4. **§8.3 中 $G$ 的定义式，讲稿把 $\gamma$ 的方程写成了 $\Delta_y\gamma=0$（无负号）**，而 [[07-基本解与广义函数]] 里 Γ 满足的是 $-\Delta\Gamma=\delta_0$。两者不矛盾（$\gamma$ 只要求调和，符号无影响），但同一个式子在同一时期内时而写 $-\Delta$、时而写 $\Delta$，容易在代入 Green 公式时出错。本篇统一为：$\gamma$ 调和（$\Delta_y\gamma=0$），$-\Delta\Gamma=\delta_0$。

---

系列：[[07-基本解与广义函数]] · [[09-Green函数的性质]] · [[10-球与半空间上的Green函数]]
