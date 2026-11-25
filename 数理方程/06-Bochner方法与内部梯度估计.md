# 06 Bochner 方法与内部梯度估计

> 系列：[[04-调和函数·Harnack不等式与内部梯度估计]] · [[05-下调和函数·弱极值原理与Hopf边界点引理]] · [[07-基本解与广义函数]]

[[04-调和函数·Harnack不等式与内部梯度估计]] 给出的估计

$$
|\nabla u(0)|\le\frac{n}{R}u(0)
$$

有显式常数，代价是必须假设 $u\ge0$。本篇的 **Bochner 方法** 是一条互补的路线：**不要求 $u\ge0$**，代价是常数只知道依赖维数。它的工具箱是「算 $\Delta|\nabla u|^2$ + 截断函数 + 弱极值原理」，这是椭圆方程先验估计的标准套路。

---

## 6.1 定理与证明

**定理（内部梯度估计，Bochner 方法）**　设 $u\in C^3(B_1)\cap C(\overline{B_1})$ 调和，则

$$
\sup_{B_{1/2}}|\nabla u|\le C_n\sup_{\partial B_1}|u|,
$$

其中 $C_n$ 只依赖于维数 $n$。

**证明**　要估计的量是 $|\nabla u|^2$，所以先算它的 Laplacian。

### 第一步：Bochner 恒等式

$$
\Delta\bigl(|\nabla u|^2\bigr)
=\sum_{i=1}^n\bigl(|\nabla u|^2\bigr)_{x_ix_i}
=\sum_{i=1}^n\left(\sum_{k=1}^n(u_{x_k})^2\right)_{x_ix_i}
=\sum_{i,k=1}^n2\bigl(u_{x_kx_i}u_{x_k}\bigr)_{x_i}
=2\sum_{i,k=1}^n\Bigl[u_{x_ix_ix_k}u_{x_k}+(u_{x_ix_k})^2\Bigr]
=2\langle\nabla(\Delta u),\nabla u\rangle+2|\nabla^2u|^2 .
$$

（其中 $|\nabla^2u|^2=\displaystyle\sum_{i,k}(u_{x_ix_k})^2$ 为 Hessian 的 Frobenius 范数的平方。第三行是先对 $\partial_{x_i}$ 求一次导、把 $\sum_k2u_{x_k}u_{x_kx_i}$ 写成 $\sum_k2(u_{x_kx_i}u_{x_k})$ 再求一次导；第四行即 Leibniz 展开。）

**代入方程 $\Delta u=0$**，得**关键不等式**

$$
\Delta\bigl(|\nabla u|^2\bigr)=2|\nabla^2u|^2\ge0 .
$$

即 $|\nabla u|^2$ **是下调和函数**。但它本身在 $B_1$ 上没有界（$\partial B_1$ 上可能很大），所以不能直接对它用极值原理——必须乘一个在边界附近消失的截断函数。

### 第二步：构造截断函数

令 $\varphi$ 为**截断函数**（cutoff），$\varphi\in C_0^\infty(B_1)$，满足

$$
\varphi\equiv1\ \text{in }B_{1/2},\qquad
\varphi\ge0\ \text{in }B_{3/4}\setminus\overline{B_{1/2}},\qquad
\varphi\equiv0\ \text{in }B_1\setminus B_{3/4},
$$

并且

$$
0\le\varphi\le1,\qquad|\nabla\varphi|^2+|\Delta\varphi|\le c_n,
$$

其中 $c_n$ 只依赖于维数。这样的函数存在，可以利用 $e^{-1/t}$ 构造（见 [[12-习题课]] 第二次作业第 1 题）。

### 第三步：算 $\Delta(|\nabla u|^2\varphi^2)$

$$
\Delta\bigl(|\nabla u|^2\varphi^2\bigr)
=\sum_{i=1}^n\bigl(|\nabla u|^2\varphi^2\bigr)_{x_ix_i}
=\sum_{i=1}^n\Bigl[(|\nabla u|^2)_{x_i}\varphi^2+|\nabla u|^2\cdot2\varphi\varphi_{x_i}\Bigr]_{x_i}
=\Delta\bigl(|\nabla u|^2\bigr)\varphi^2+4\varphi\bigl\langle\nabla\bigl(|\nabla u|^2\bigr),\nabla\varphi\bigr\rangle
+2|\nabla u|^2|\nabla\varphi|^2+2\varphi|\nabla u|^2\Delta\varphi .
$$

由 $\nabla(|\nabla u|^2)=2\nabla^2u\cdot\nabla u$ 以及 **Young 不等式** $2ab\le a^2+b^2$（取 $a=\varphi|\nabla^2u|$、$b=4|\nabla u||\nabla\varphi|$）：

$$
4\varphi\bigl\langle\nabla\bigl(|\nabla u|^2\bigr),\nabla\varphi\bigr\rangle
=8\varphi\bigl\langle\nabla^2u\,\nabla u,\nabla\varphi\bigr\rangle
\ge-8\varphi\,|\nabla^2u|\,|\nabla u|\,|\nabla\varphi|
\ge-\varphi^2|\nabla^2u|^2-16|\nabla u|^2|\nabla\varphi|^2 .
$$

因此

$$
|\nabla^2u|^2\varphi^2+4\varphi\bigl\langle\nabla\bigl(|\nabla u|^2\bigr),\nabla\varphi\bigr\rangle+16|\nabla u|^2|\nabla\varphi|^2\ge0 .
$$

代入第三步的表达式，并利用第一步的 $2|\nabla^2u|^2$：

$$
\Delta\bigl(|\nabla u|^2\varphi^2\bigr)
\ge2|\nabla^2u|^2\varphi^2-\varphi^2|\nabla^2u|^2-16|\nabla u|^2|\nabla\varphi|^2+2|\nabla u|^2|\nabla\varphi|^2+2\varphi|\nabla u|^2\Delta\varphi
\ge\varphi^2|\nabla^2u|^2-\bigl(14|\nabla\varphi|^2+2|\varphi\Delta\varphi|\bigr)|\nabla u|^2
\ge-C_n|\nabla u|^2\qquad\text{于 }B_1 .
$$

### 第四步：配上一个「补偿项」并应用弱极值原理

注意

$$
\Delta(u^2)=2u\Delta u+2|\nabla u|^2=2|\nabla u|^2\qquad(\text{用 }\Delta u=0).
$$

于是

$$
\Delta\bigl(|\nabla u|^2\varphi^2+C_nu^2\bigr)\ge-C_n|\nabla u|^2+C_n\cdot2|\nabla u|^2=C_n|\nabla u|^2\ge0\qquad\text{于 }B_1 .
$$

**这是一个下调和函数。** 又由于 $\varphi$ 在 $\partial B_1$ 的邻域内恒为零，函数

$$
F:=|\nabla u|^2\varphi^2+C_nu^2
$$

可以连续延拓到 $\overline{B_1}$，且在 $\partial B_1$ 上 $F=C_nu^2$。由 [[05-下调和函数·弱极值原理与Hopf边界点引理]] 的弱极值原理（对一般有界区域的形式），

$$
\sup_{B_1}\bigl(|\nabla u|^2\varphi^2+C_nu^2\bigr)
\le\sup_{\partial B_1}\bigl(|\nabla u|^2\varphi^2+C_nu^2\bigr)
=C_n\sup_{\partial B_1}u^2 .
$$

**在 $B_{1/2}$ 上 $\varphi\equiv1$**，故对 $x\in B_{1/2}$，

$$
|\nabla u(x)|^2\le|\nabla u(x)|^2+C_nu(x)^2\le C_n\sup_{\partial B_1}u^2 ,
$$

取上确界并开平方即得

$$
\sup_{B_{1/2}}|\nabla u|\le\sqrt{C_n}\sup_{\partial B_1}|u| .
$$

$\square$

> **补充（为什么用 $\varphi^2$ 而不是 $\varphi$）**　这是全篇最值得琢磨的一处。主项来自 $\Delta(|\nabla u|^2)$，即 $2|\nabla^2u|^2\varphi^2$；交叉项 $4\varphi\langle\nabla(|\nabla u|^2),\nabla\varphi\rangle$ 用 Young 不等式之后会产生**一项 $-\varphi^2|\nabla^2u|^2$**。注意它恰好是 $\varphi^2$（而不是 $\varphi$）——因为 Young 的 $a$ 取了 $\varphi|\nabla^2u|$、$b$ 取了 $4|\nabla u||\nabla\varphi|$，其中 $b$ 自带一个 $\varphi$ 的幂被吸收。于是 $2|\nabla^2u|^2\varphi^2-\varphi^2|\nabla^2u|^2=|\nabla^2u|^2\varphi^2\ge0$，**主项被完整地「吃掉」了**，只剩下 $|\nabla u|^2$ 的项。

$$
\text{若改用 }\varphi:\quad
\Delta\bigl(|\nabla u|^2\varphi\bigr)\ge|\nabla^2u|^2\varphi-\bigl(\cdots\bigr)|\nabla u|^2,
$$

> 主项变成 $\varphi|\nabla^2u|^2$ 与 $|\nabla^2u|^2\varphi$ 相减——**还能配上，但交叉项里会多出 $\varphi^{-1}$ 型系数**，因为要写成 $\varphi\cdot(\varphi|\nabla^2u|^2)$ 时 $b^2$ 里的 $|\nabla\varphi|^2$ 无法被 $\varphi$ 控制（$|\nabla\varphi|/\varphi$ 在截断的过渡区无界）。所以 $\varphi^2$ 是让「主项配平」与「交叉项可被截断函数的导数控制」同时成立的最小幂次。

> **补充（Bochner 恒等式的地位）**　

$$
\Delta|\nabla u|^2=2|\nabla^2u|^2+2\langle\nabla(\Delta u),\nabla u\rangle
$$

> 对**任何** $C^3$ 函数成立，称为 **Bochner 恒等式**。当 $\Delta u=0$ 时右端第一项非负——这就是「调和函数的梯度平方是下调和函数」这一事实的来源，也是后续一切梯度估计的起点。对一般的 $Lu=\sum a_{ij}u_{x_ix_j}+\sum b_iu_{x_i}=0$（$a_{ij}$ 正定，见 [[12-习题课]] 第二次作业第 7 题）也有一条加权版本的 Bochner 恒等式，做法是把 $|\nabla u|^2$ 换成关于 $a_{ij}$ 的二次型 $\sum a_{ij}u_{x_i}u_{x_j}$。

> **补充（与 [[04-调和函数·Harnack不等式与内部梯度估计]] 的分工）**　

| 路线             | 结论                                                               | 是否需要 $u\ge0$ | 常数  |
| -------------- | ---------------------------------------------------------------- | --- | -------- |
| 方向导数 + 散度定理    | $\nabla u(0)\le\dfrac{n}{R}u(0)$                                 | 需要  | 显式 $n/R$ |
| Bochner + 截断函数 | $\sup_{B_{1/2}}\nabla u\le C_n\sup_{\partial B_1}\lvert u\rvert$ | 不需要 | 只知依赖 $n$ |

> 前者用于 Liouville 定理（要把 $R\to\infty$，必须常数显式）；后者用于缺少下界信息的场合，例如证明 $|u(x)|=o(|x|^{k+1})\Rightarrow u$ 为多项式（[[12-习题课]] 第二次作业第 8 题）。

---

## 订正

1. **定理与证明中常数 $C_n$ 的实际值不同。** 定理写的是 $\sup_{B_{1/2}}|\nabla u|\le C_n\sup_{\partial B_1}|u|$，而证明走的是先估 $|∇u|²$、最后开平方的路线，实际得到的是 $\sup|\nabla u|\le\sqrt{C_n}\sup|u|$。由于两处都只声明「只依赖 $n$」，结论不受影响；但这是一个同一符号在不同位置取不同值的复用，引用常数时需注意已经开过平方。
2. **第三步中 $|\nabla u|^2\varphi^2$ 关于 $x_i$ 求二阶导的展开，讲稿未标出括号位置。** 讲稿该行印作 $\Delta(|\nabla u|^2)=\sum_{i,k}2(u_{x_kx_i},u_{x_k})_{x_i}$（经核对原页与源文件，括号是印出来的、只是纯文本抽取时会丢失）。此式**是对的**：它把 $\sum_k2u_{x_k}u_{x_kx_i}$ 写成 $\sum_k2(u_{x_kx_i}u_{x_k})$ 后再对 $x_i$ 求一次导，与下一行 $2\sum[u_{x_ix_ix_k}u_{x_k}+(u_{x_ix_k})^2]$ 一致。本篇保留原式并补了括号，**此处不是错误**。
3. **「$|\nabla u|^2\varphi^2$ 在 $B_1$ 内取到最大值（思考为何用 $\varphi^2$？）」这句话的作用容易被误读。** 证明并不真的对 $|\nabla u|^2\varphi^2$ 用极值原理（它在 $B_1$ 上是否取到最大值无关紧要），而是对 $|\nabla u|^2\varphi^2+C_nu^2$ 用弱极值原理。括号里的「为何用 $\varphi^2$」是针对**配平主项**的，答案见本篇补充。
4. **定理的光滑性假设 $\Delta(|\nabla u|^2\varphi^2)$ 可算需要 $u\in C^3$**（$|\nabla u|^2$ 的二阶导 = $u$ 的三阶导），而非 $C^2$。讲稿写的是 $u\in C^3(B_1)\cap C(\overline{B_1})$，正确，此处只是提醒不要把 $C^3$ 误当成多余假设。

---

系列：[[04-调和函数·Harnack不等式与内部梯度估计]] · [[05-下调和函数·弱极值原理与Hopf边界点引理]] · [[07-基本解与广义函数]]
