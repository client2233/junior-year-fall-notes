# 04 调和函数·Harnack 不等式与内部梯度估计

> 系列：[[03-调和函数·平均值性质与极值原理]] · [[05-下调和函数·弱极值原理与Hopf边界点引理]] · [[06-Bochner方法与内部梯度估计]]

本篇给出调和函数的两个**定量**性质。它们都从平均值性质出发，但结论比极值原理强得多：

- **Harnack 不等式**：非负调和函数在紧子区域上的**相对**波动被常数控制；
- **内部梯度估计**：非负调和函数的一阶导数被函数值本身控制，进而给出 Liouville 定理。

---

## 4.1 Harnack 不等式

**定理（Harnack 不等式）**　令 $\Omega\subseteq\mathbb R^n$ 为一区域。若 $u\in C^2(\Omega)$ 为非负调和函数，则对任意子区域 $\Omega'\Subset\Omega$，存在常数 $C(\Omega',\Omega)>0$，使得

$$
\frac{1}{C}u(y)\le u(x)\le Cu(y),\qquad\forall x,y\in\overline{\Omega'} .
$$

**证明**　当 $\Omega\ne\mathbb R^n$ 时，令

$$
R=\frac18\operatorname{dist}(\Omega',\partial\Omega)>0;
$$

当 $\Omega=\mathbb R^n$ 时，取 $R=1$。对任意 $x\in\Omega'$，先在 $B_R(x)$ 中证明。

**第一步：半径 $R$ 内的比较。** 任取 $y\in B_R(x)$，由于

$$
B_R(x)\subset B_{2R}(y)\subset B_{3R}(x)\Subset\Omega
$$

（前一个包含用 $|y-x|<R$，后一个用 $|y-x|<R$，最后一个用 $3R<\operatorname{dist}(\Omega',\partial\Omega)$），由平均值性质与 $u\ge0$（扩大积分区域不会使积分变小）得

$$
u(y)=\frac{1}{|B_{2R}(y)|}\int_{B_{2R}(y)}u(w)\,dw
=\frac{1}{2^n|B_R(x)|}\int_{B_{2R}(y)}u(w)\,dw
\ge\frac{1}{2^n|B_R(x)|}\int_{B_R(x)}u(w)\,dw
=\frac{1}{2^n}u(x).
$$

任取 $z\in B_R(x)$，同理有

$$
u(x)\ge\frac{1}{2^n}u(z).
$$

因此

$$
u(y)\ge\frac{1}{4^n}u(z),\qquad\forall y,z\in B_R(x).
$$

**第二步：用球链拼到整个 $\Omega'$。** 由 $\overline{\Omega'}$ 的紧性，可取 $x_1,\dots,x_N\in\Omega'$，使得 $B_i=B_R(x_i)$ 覆盖 $\overline{\Omega'}$。因 $\Omega'$ 连通且每个 $B_i$ 均与 $\Omega'$ 相交，这些球的并集连通；任意两点可由至多 $N$ 个相邻相交的球连接，且这些球均包含于 $\Omega$。沿球链逐次比较（每跨一步损失因子 $4^n$），即可取

$$
C(\Omega',\Omega)=4^{nN}=2^{2nN},
$$

其中 $N$ 只依赖于 $\Omega'$ 与 $\Omega$。$\square$

**直观理解**　对

$$
\begin{cases}
-\Delta u=0, & \text{in }\Omega,\\
u=g, & \text{on }\partial\Omega
\end{cases}
$$

来说，$g$ 的波动可以非常大、任意大，但只要**往里走**，$u$ 的相对波动就会被控制住。

> **要点**　Harnack 不等式与极值原理的差别在于「绝对」与「相对」：极值原理说的是 $u\le\max_{\partial\Omega}g$（绝对上界由边界数据给出，可以很大，而且依赖具体边界值）；Harnack 说的是 $u(x)/u(y)\le C$（管内两点的比值被一个**只依赖几何**的常数控制，与边界数据的尺度完全无关）。后者的实用性正是「先验估计」的雏形。

> **补充（常数为什么只能依赖 $\Omega',\Omega$）**　球链的步数 $N$ 由 $\overline{\Omega'}$ 的覆盖数决定，每步的损失 $4^n$ 由维数决定——两者都不含 $u$ 或 $g$。这就是 Harnack 常数的「普适性」。反过来，当 $x,y$ 靠近 $\partial\Omega$ 时，$R\to0$、步数 $N\to\infty$，常数就会发散，这与边界处 $u$ 可以任意剧烈波动是一致的。

---

## 4.2 内部梯度估计

**定理（内部梯度估计）**　记 $B_R=B_R(0)$。设 $u\in C^3(\overline{B_R})\cap C^1(\overline{B_R})$ 为**非负**调和函数，则

$$
|\nabla u(0)|\le\frac{n}{R}\,u(0).
$$

**证明**　关键观察：由于 $\Delta u=0$ 且 $\Delta$ 是常系数算子，

$$
\Delta\frac{\partial u}{\partial x_i}=0,\qquad i=1,2,\dots,n,
$$

即**调和函数的任意一阶偏导数仍调和**。

**第一步：用旋转不变性把梯度化成一个方向导数。** 令 $y=Ax$，其中 $A=(A_{ij})$ 为正交矩阵（$A^TA=I$），并令 $\tilde u(y)=u(A^{-1}y)$。由于 $A^{-1}=A^T$，有

$$
\frac{\partial\tilde u}{\partial y_i}(y)=\sum_{j=1}^n\frac{\partial x_j}{\partial y_i}\frac{\partial u}{\partial x_j}(x)=\sum_{j=1}^nA_{ij}u_{x_j}(x),\qquad x=A^{-1}y .
$$

再求一次导：

$$
\frac{\partial^2\tilde u}{\partial y_i^2}(y)=\sum_{j,k=1}^nA_{ij}A_{ik}u_{x_jx_k}(x),
$$

于是

$$
\Delta_y\tilde u(y)=\sum_{i,j,k=1}^nA_{ij}A_{ik}u_{x_jx_k}(x)
=\sum_{j,k=1}^n\delta_{jk}u_{x_jx_k}(x)
=\sum_{j=1}^nu_{x_jx_j}(x)=(\Delta u)(Ax)=0,
$$

其中用了正交性 $\displaystyle\sum_{i=1}^nA_{ij}A_{ik}=\delta_{jk}$。**Laplace 算子在旋转下不变。**

适当旋转坐标（并把旋转后的函数仍记为 $u$），使

$$
|\nabla u(0)|=\frac{\partial u}{\partial x_1}(0)\ \ge0.
$$

**第二步：把方向导数写成平均值。** 由第一步，$\dfrac{\partial u}{\partial x_1}$ 仍调和，故由平均值性质（体积形式）

$$
\frac{\partial u}{\partial x_1}(0)=\frac{1}{|B_R|}\int_{B_R}\frac{\partial u}{\partial x_1}\,dx .
$$

**第三步：用散度定理把体积分化成边界通量。** 考虑向量值函数

$$
X=(u,0,\dots,0),
$$

则

$$
\operatorname{div}X=\frac{\partial u}{\partial x_1}+0+\cdots+0,
$$

从而

$$
\int_{B_R}\frac{\partial u}{\partial x_1}\,dx=\int_{B_R}\operatorname{div}X\,dx=\int_{\partial B_R}\langle X,\nu\rangle\,d\sigma_x=\int_{\partial B_R}u\,\nu_1\,d\sigma_x .
$$

**第四步：用非负性与 $|\nu_1|\le1$ 收尾。**

$$
\frac{\partial u}{\partial x_1}(0)=\frac{1}{|B_R|}\int_{\partial B_R}u\,\nu_1\,d\sigma_x
\le\frac{1}{|B_R|}\int_{\partial B_R}u\,|\nu_1|\,d\sigma_x
\le\frac{1}{|B_R|}\int_{\partial B_R}u\,d\sigma_x
=\frac{|\partial B_R|}{|B_R|}u(0)
=\frac{n}{R}u(0),
$$

最后一步用了 $|\partial B_R|/|B_R|=n/R$。结合 $|\nabla u(0)|=\partial_{x_1}u(0)$ 即得结论。$\square$

### 4.2.1 三条观察

**观察 1（$0$ 点不重要）**　对非负调和函数 $u\in C^3(\overline{B_R(x_0)})\cap C^1(\overline{B_R(x_0)})$，令 $\hat u(x):=u(x+x_0)$，则

$$
|\nabla\hat u(0)|=|\nabla u(x_0)|\le\frac{n}{R}u(x_0).
$$

更一般地，设 $u\in C^3(\Omega)$ 为非负调和函数，$\Omega'\Subset\Omega$，且 $0<R<\operatorname{dist}(\Omega',\partial\Omega)$，则

$$
|\nabla u(x)|\le\frac{n}{R}u(x),\qquad x\in\Omega' .
$$

**观察 2（Laplace 算子旋转不变）**　即证明第一步中所算的：对 $A\in\mathbb R^{n\times n}$、$A^TA=I$ 与 $\hat u(x):=u(Ax)$，有

$$
\Delta\hat u(x)=\sum_{j,k=1}^n\delta_{jk}\frac{\partial^2u}{\partial x_j\partial x_k}(Ax)=(\Delta u)(Ax).
$$

**观察 3（推论：Liouville 定理）**　若 $u\in C^3(\mathbb R^n)$ 为调和函数，且存在常数 $C$ 使得 $u\ge-C$，则 $u$ 为常数。

**证明**　令 $\hat u(x)=u(x)+C\ge0$，则 $\hat u$ 在 $\mathbb R^n$ 上非负调和。对任意 $x_0\in\mathbb R^n$，

$$
|\nabla u(x_0)|=|\nabla\hat u(x_0)|\le\frac{n}{R}\bigl(u(x_0)+C\bigr)\xrightarrow[\ R\to+\infty\ ]{}0,
$$

故 $\nabla u(x_0)=0$。由 $x_0$ 任意，$u$ 为常数。$\square$

> **补充（非负性不可省）**　$u(x)=x_1$ 在 $\mathbb R^n$ 上调和，$|\nabla u|\equiv1$ 不随 $R$ 衰减。它之所以不违背梯度估计，正是因为它**没有单侧下界**。Liouville 定理的假设也只能是「单侧有界」而不是「有界」——单侧即可，因为 $\hat u=u+C\ge0$ 就把单侧界用掉了。

> **补充（两种梯度估计的比较）**　本篇的估计是 $|\nabla u(0)|\le\dfrac{n}{R}u(0)$，两个特点：常数 $n/R$ 显式、但要求 $u\ge0$，且只控制**球心**处的梯度。[[06-Bochner方法与内部梯度估计]] 用 Bochner 技巧给出

$$
\sup_{B_{1/2}}|\nabla u|\le C_n\sup_{\partial B_1}|u|,
$$

它**不需要** $u\ge0$，代价是常数不显式（只依赖维数），且用到了截断函数与弱极值原理。两条路线各有用途：要显式常数（如证明 Liouville）用前者，要摆脱非负假设用后者。

> **补充（从梯度估计还能得到什么）**　反复对 $\partial_{i_1}\cdots\partial_{i_{j-1}}u$（它也调和）用一阶估计，可把 $k$ 阶导数用 $u$ 本身控制：

$$
|\partial_{i_1}\cdots\partial_{i_k}u(x_0)|\le\left(\frac{2kn}{R}\right)^k\sup_{B_R(x_0)}|u|.
$$

由此可证：若 $|u(x)|=o(|x|^{k+1})$，则 $u$ 是次数不超过 $k$ 的多项式（**调和多项式**）。完整推导见 [[12-习题课]] 第二次作业第 8 题。

---

## 订正

1. **Harnack 常数的形式存疑。** 讲稿印的是 $C(\Omega',\Omega)=2^{2N}n$。而按证明里的球链论证，每跨一个球给出因子 $4^n=2^{2n}$，$N$ 个球给出 $4^{nN}=2^{2nN}$——**$n$ 应当出现在指数上**。两种写法都只依赖 $\Omega'$ 与 $\Omega$（不影响不等式成立与「常数只依赖几何」这一结论），但若需要写出常数的具体取值，应按 $2^{2nN}$ 理解。本篇正文采用 $4^{nN}=2^{2nN}$。
2. **Harnack 定理的适用范围写成 $\overline{\Omega'}$ 更准确。** 讲稿正文与证明处时而写 $\Omega'$、时而写 $\overline{\Omega'}$；证明第二步明确用了「$\overline{\Omega'}$ 的紧性」，故不等式的成立范围应是 $\overline{\Omega'}$（讲稿在 $\Omega'$ 上加的一横在抽取时容易丢失）。本篇已统一为 $\overline{\Omega'}$。
3. **内部梯度估计的光滑性假设。** 讲稿（两处）均写 $u\in C^3(B_R)\cap C^1(\overline{B_R})$。证明中「$\partial u/\partial x_1$ 仍调和」进而「可用平均值性质」需要 $\partial u/\partial x_1\in C^2$，即 $u\in C^3$（本轮只用在球内部，边界上只需 $C^1$ 保证 $u$ 在 $\partial B_R$ 上取值）。本篇写作 $u\in C^3(\overline{B_R})\cap C^1(\overline{B_R})$，与讲稿等价但避免了记号歧义。
4. **第三步的不等号方向。** 讲稿直接写 $\displaystyle\int_{\partial B_R}u\nu_1d\sigma\le\int_{\partial B_R}u|\nu_1|d\sigma$ 是成立的（因 $u\ge0$），但要点在于 $\partial_{x_1}u(0)$ 可能是负的——必须**先旋转坐标使 $|\nabla u(0)|=\partial_{x_1}u(0)\ge0$**，这一步「可做到」的理由（取 $e_1$ 沿 $\nabla u(0)$ 方向）讲稿未展开。本篇已在正文补上。

---

系列：[[03-调和函数·平均值性质与极值原理]] · [[05-下调和函数·弱极值原理与Hopf边界点引理]] · [[06-Bochner方法与内部梯度估计]]
