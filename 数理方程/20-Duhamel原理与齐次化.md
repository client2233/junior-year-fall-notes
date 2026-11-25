# 20 Duhamel 原理（齐次化原理）

> 系列：[[19-分离变量解的正则性·三个收敛命题]] · [[18-热方程·建模、边值条件与分离变量法]] · [[21-Fourier变换与速降函数空间]]

[[19-分离变量解的正则性·三个收敛命题]] 解决了**齐次方程 + 非零初值**。现实中还有另一类问题：**方程本身带源项**（体内有热源）。这一篇讲怎样把这种问题化归成已经会解的那种。

结论一句话：**带源项的问题 = 齐次问题的「初值」沿时间连续叠加。**

---

## 20.1 带热源的建模

物理场景：$t>0$ 时刻在 $x$ 处注入热量使其升温（或降温），记热源强度为 $f(x,t)$。重做一遍 [[18-热方程·建模、边值条件与分离变量法]] §18.1 的能量收支。

任取 $V\Subset\Omega$，取 $t_1<t_2$。$V$ 内能量的增量来自两处：**穿过边界流入的热量**与**体内热源产生的热量**。写成「$t_1$ 时刻能量减 $t_2$ 时刻能量」（差一个负号，是为了让下面两式的符号自然对上）：

$$
E_V(t_1)-E_V(t_2)=\int_{t_2}^{t_1}\!\!\int_{\partial V}\langle\vec q(x,t),-\nu\rangle\,d\sigma_x\,dt+\int_{t_2}^{t_1}\!\!\int_Vf(x,t)\,dx\,dt .
$$

另一方面，由能量与温度的关系 $E_V(t)=c\int_Vu(x,t)dx$（这里 $c$ 已把 $\rho$ 吸收进去，见 [[18-热方程·建模、边值条件与分离变量法]] 的订正第 3 条），

$$
E_V(t_1)-E_V(t_2)=\int_{t_2}^{t_1}\left[\frac{d}{dt}\int_Vc\,u(x,t)\,dx\right]dt .
$$

两式相等。注意 $\int_{t_2}^{t_1}=-\int_{t_1}^{t_2}$，两边同乘 $-1$ 并把 $\vec q=-k\nabla u$ 代入：

$$
\int_{t_1}^{t_2}\!\!\int_Vc\,u_t\,dx\,dt=\int_{t_1}^{t_2}\!\!\int_V\left(k\Delta u+f\right)dx\,dt .
$$

$V$ 与 $[t_1,t_2]$ 任取，故被积函数相等。取 $c=k=1$（与 [[18-热方程·建模、边值条件与分离变量法]] 同一归一）：

$$
\boxed{\ u_t-\Delta u=f(x,t)\ }\tag{20.1}
$$

连同初边值，得到完整的**非齐次初边值问题**：

$$
\begin{cases}u_t-\Delta u=f(x,t), & (x,t)\in\Omega\times(0,+\infty),\\[4pt] u(x,0)=\varphi(x), & x\in\Omega,\\[4pt] u(x,t)=0, & x\in\partial\Omega,\ t>0 .\end{cases}\tag{20.2}
$$

> **要点**　$f$ 与 $\varphi$ 的地位完全不同：$\varphi$ 只作用一次（$t=0$），$f$ 在**每个时刻**都在作用。所以解法上要把 $f$ 沿时间「铺开」成一族初值 —— 这正是 Duhamel 原理。

---

## 20.2 原理的陈述

问题 (20.2) 的解等于下面两个问题的解之和。

**方程 1（齐次方程，非零初值）**：

$$
\begin{cases}\partial_tu_1-\Delta u_1=0, & (x,t)\in\Omega\times(0,+\infty),\\[4pt] u_1(x,0)=\varphi(x), & x\in\Omega,\\[4pt] u_1(x,t)=0, & (x,t)\in\partial\Omega\times(0,+\infty).\end{cases}
$$

**方程 2（对每个「注入时刻」$\tau>0$ 一个齐次问题）**：对每个 $\tau>0$，记 $w_\tau(x,s)$ 为

$$
\begin{cases}\partial_sw_\tau-\Delta w_\tau=0, & (x,s)\in\Omega\times(0,+\infty),\\[4pt] w_\tau(x,0)=f(x,\tau), & x\in\Omega,\\[4pt] w_\tau(x,s)=0, & (x,s)\in\partial\Omega\times(0,+\infty)\end{cases}\tag{20.3}
$$

的解，其中 $s=t-\tau$ 是「从注入时刻起算的时间」。

**则 (20.2) 的解为**

$$
u(x,t)=u_1(x,t)+\int_0^t w_\tau(x,t-\tau)\,d\tau .\tag{20.4}
$$

> **要点（为什么 $f$ 变成了「初值」）**　热方程对时间是**一阶**的（[[18-热方程·建模、边值条件与分离变量法]] §18.1 的对照表）。一阶方程只有「初值」，没有「右端项延拓」这回事；而源项 $f(x,\tau)$ 在 $\tau$ 时刻注入的那一份热量，从 $\tau$ 之后就跟「初始温度」一样按齐次方程自由扩散。所以把每一份注入都当成一个**以 $f(\cdot,\tau)$ 为初值、从 $\tau$ 起跑**的齐次问题，(20.4) 就是把它们的贡献积分起来。
>
> 反过来：**如果方程对时间是二阶的（波方程），同一套话仍成立，但换元后要积分 $\int_0^t\sin(\cdots)$**，见 §20.5。

---

## 20.3 验证

只需验证 $u_2(x,t):=\displaystyle\int_0^tw_\tau(x,t-\tau)d\tau$ 满足

$$
\begin{cases}\partial_tu_2-\Delta u_2=f(x,t), & (x,t)\in\Omega\times(0,+\infty),\\[4pt] u_2(x,0)=0, & x\in\Omega,\\[4pt] u_2(x,t)=0, & (x,t)\in\partial\Omega\times(0,+\infty).\end{cases}\tag{20.5}
$$

（记 (20.5) 为**方程 3**。）然后由线性性，$u_1+u_2$ 满足 (20.2)。

**证明（$\Leftarrow$）**　先算时间导数。由 Leibniz 法则（上限依赖 $t$、被积函数也依赖 $t$）：

$$
\partial_tu_2(x,t)=w_t(x,t-t)\cdot1+\int_0^t\partial_tw_\tau(x,t-\tau)\,d\tau .
$$

第一项是上限项，$\tau=t$ 时 $t-\tau=0$，故等于 $w_t(x,0)=f(x,t)$。于是

$$
\begin{aligned}
\left(\partial_t-\Delta\right)u_2(x,t)
&=w_t(x,0)+\int_0^t\partial_tw_\tau(x,t-\tau)\,d\tau-\int_0^t\Delta w_\tau(x,t-\tau)\,d\tau\\
&=w_t(x,0)+\int_0^t\underbrace{\left[\partial_sw_\tau-\Delta w_\tau\right]}_{=0\ \text{由 (20.3)}}(x,t-\tau)\,d\tau\\
&=w_t(x,0)=f(x,t).
\end{aligned}
$$

（第二行里把「对 $t$ 求导」换成「对 $s=t-\tau$ 求导」，因为 $\partial_t w_\tau(x,t-\tau)=\partial_sw_\tau(x,t-\tau)$。）

初始条件：$u_2(x,0)=\int_0^0(\cdots)=0$。边值：每个 $w_\tau(\cdot,s)$ 在 $\partial\Omega$ 上为零，故它们的积分也为零。于是 $u_2$ 解 (20.5)，$u_1+u_2$ 解 (20.2)。$\square$

**证明（$\Rightarrow$）**　设 $u$ 解 (20.2)，则 $u-u_1$ 满足 (20.5)（齐次边值、零初值、右端 $f$）。所以 (20.4) 不是「造出了一个解」，而是**唯一解**（唯一性由 [[18-热方程·建模、边值条件与分离变量法]] §18.5 的能量估计给出）。$\square$

> **要点（上限项就是源项本身）**　整个证明只有一处实质：Leibniz 法则给出的那一个**上限项** $w_t(x,0)=f(x,t)$。它之所以等于 $f(x,t)$ 而不是别的，是因为 (20.3) 的初值就是 $f(\cdot,\tau)$，而在 $\tau=t$ 这个「最新的注入」处 $s=t-\tau=0$。**如果源项不在初值里而在方程右端，这一项就出不来。**

---

## 20.4 一维情形的显式解

考虑

$$
\begin{cases}u_t-u_{xx}=f(x,t), & (x,t)\in(0,L)\times(0,+\infty),\\[4pt] u(x,0)=0, & x\in(0,L),\\[4pt] u(0,t)=u(L,t)=0, & t>0 .\end{cases}
$$

对每个 $\tau>0$，解 (20.3)。原问题的 $C^0$ 相容性条件变为

$$
f(0,\tau)=f(L,\tau)=0\qquad(\tau>0).
$$

否则 $w_\tau$ 在角点 $(0,0)$ 处就不连续（同 [[19-分离变量解的正则性·三个收敛命题]] §19.1）。

**解 $w_\tau$。** 记 $f(\cdot,\tau)$ 的正弦展开为

$$
f(x,\tau)=\sum_{k=1}^{\infty}a_k(\tau)\sin\frac{k\pi x}{L},\qquad a_k(\tau)=\frac2L\int_0^Lf(y,\tau)\sin\frac{k\pi y}{L}\,dy .
$$

由 [[18-热方程·建模、边值条件与分离变量法]] §18.4，把初值 $f(\cdot,\tau)$ 换成上面的展开：

$$
w_\tau(x,s)=\sum_{k=1}^{\infty}a_k(\tau)e^{-\lambda_ks}\sin\frac{k\pi x}{L},\qquad s=t-\tau>0 .
$$

**叠加。** 代入 (20.4)：

$$
u(x,t)=\int_0^tw_\tau(x,t-\tau)\,d\tau=\sum_{k=1}^{\infty}\left[\int_0^ta_k(\tau)e^{-\lambda_k(t-\tau)}\,d\tau\right]\sin\frac{k\pi x}{L}.
$$

再把 $a_k(\tau)$ 展开回去：

$$
u(x,t)=\sum_{k=1}^{\infty}\frac2L\int_0^L\left[\int_0^tf(y,\tau)e^{-\lambda_k(t-\tau)}\,d\tau\right]\sin\frac{k\pi y}{L}\,dy\ \sin\frac{k\pi x}{L}.
$$

这就是一维非齐次问题的解。

> **要点（读这个公式的方法）**　括号里的 $\int_0^ta_k(\tau)e^{-\lambda_k(t-\tau)}d\tau$ 是**一个卷积**：$a_k$ 在时刻 $\tau$ 注入，之后按 $e^{-\lambda_k(t-\tau)}$ 衰减到 $t$ 时刻。$\lambda_k$ 越大衰减越快，所以**高频的源项「来即走」，低频的源项会积累**。若 $f$ 与 $t$ 无关，则积分给出 $\dfrac{a_k}{\lambda_k}\left(1-e^{-\lambda_kt}\right)\to\dfrac{a_k}{\lambda_k}$（$t\to\infty$），即稳态解 $u_*=\sum_k\frac{a_k}{\lambda_k}\sin\frac{k\pi x}{L}$，它满足 $-u_*''=f$ —— 与「长时间后温度不再变化」的物理直觉一致。

---

## 20.5 为什么叫「齐次化」：半群语言

**补充（抽象形式）**　把 $A:=\Delta$（配 Dirichlet 边值条件）看成一个算子，$e^{tA}\varphi$ 就是齐次问题的解算子（**热半群**）。则 (20.4) 写成

$$
u(t)=e^{tA}\varphi+\int_0^te^{(t-\tau)A}f(\tau)\,d\tau .
$$

这是一阶线性方程 $u'=Au+g$ 的**常数变易公式**：比较常微分方程

$$
v'(t)=Av(t)+g(t),\quad v(0)=v_0\quad\Longrightarrow\quad v(t)=e^{tA}v_0+\int_0^te^{(t-\tau)A}g(\tau)\,d\tau ,
$$

两者形状完全一样。所以 Duhamel 原理不是热方程的什么特殊技巧，而是**一阶线性演化方程的通例**：$A$ 从矩阵换成 $\Delta$，矩阵指数换成热半群，公式一字不改。

**补充（这个「半群」是什么）**　$e^{tA}$ 的具体形式由特征展开给出：

$$
e^{tA}\varphi=\sum_{k=1}^{\infty}\left(\varphi,u_k\right)_{L^2}e^{-\lambda_kt}u_k ,
$$

即「在每个特征方向上乘 $e^{-\lambda_kt}$」。于是 (20.4) 又能读成：**源项在每个特征方向上的贡献，是从注入时刻起各乘一个 $e^{-\lambda_k(t-\tau)}$ 之后累积起来的**。这与 §20.4 的一维公式逐项对应。（$\lambda_k$ 与 $u_k$ 来自 [[14-谱理论·鼓面的数学建模与Laplace特征值]] 的谱定理。）

---

## 订正

1. **验证里最后一步写成 $w(x,0)=f(x,\tau)$ 是错的，应为 $f(x,t)$。** 积分 $\int_0^tw_\tau(x,t-\tau)d\tau$ 里的 $\tau$ 是**哑变量**，上限项取 $\tau=t$，此时 $t-\tau=0$，所以得到的是 $w_t(x,0)=f(x,t)$ —— 即「当前时刻的源项」。写成 $f(x,\tau)$ 是把哑变量留在了结果里，量纲上也说不通（$\tau$ 在积分外无意义）。这个错误会让读者看不懂「为什么右端恰好是 $f(x,t)$」，而这一项正是整个证明的唯一实质。
2. **「方程 2」的记号自相矛盾。** 讲稿写
$\partial_tw(x,t-\tau)-\Delta w(x,t-\tau)=0$、初值 $w(x,0)=f(x,\tau)$、边值 $w(x,t-\tau)=0$，
同一个符号 $w$ 被当成「$(x,t)$ 的函数」又当成「$(x,s)$ 的函数」，而且 $w(x,0)$ 里根本没有 $\tau$，看不出它依赖 $\tau$。本篇改用 $w_\tau(x,s)$、$s=t-\tau$，把「哪一族的解」与「族内的时间」分开。
3. **§4 最后一个等式的哑变量冲突。** 讲稿写
$u(x,t)=\sum_k\frac2L\int_0^L\left\{\int_0^tf(x,\tau)e^{-\lambda_k(t-\tau)}d\tau\right\}\sin\frac{k\pi}{L}x\,dx\ \sin\frac{k\pi}{L}x$，
内层积分（空间）与最外层的自变量都叫 $x$，外层的 $x$ 被遮蔽。应为 $\int_0^L(\cdots)dy$，内层用 $f(y,\tau)$、$\sin\frac{k\pi y}{L}$。
4. **§3 的能量式里 $c$ 的口径。** 讲稿写 $\int_Vc\,u(x,t)dx$，而 [[18-热方程·建模、边值条件与分离变量法]] 用的是 $\rho c$。两者只差一个常数，但同一门课两处不一致，容易让人以为漏了密度。本篇统一按「$c$ 已吸收 $\rho$」并取 $c=k=1$。
5. **「方程 1 的解与方程 3 的解相加，才构成原方程的解」这句话没说清为什么。** 它其实只是线性叠加：$u=u_1+(u-u_1)$，其中 $u-u_1$ 满足零初值、齐次边值、右端 $f$。讲稿把 (20.5) 单列为「方程 3」却不说明它从哪来，读者会以为这是需要额外论证的分解。本篇 §20.3 的 $(\Rightarrow)$ 部分补了这句。

---

系列：[[19-分离变量解的正则性·三个收敛命题]] · [[18-热方程·建模、边值条件与分离变量法]] · [[21-Fourier变换与速降函数空间]]
