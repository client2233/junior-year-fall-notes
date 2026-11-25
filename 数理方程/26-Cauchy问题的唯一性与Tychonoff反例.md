# 26 全空间 Cauchy 问题的唯一性 · Tychonoff 反例

> 系列：[[25-热方程的梯度估计与Harnack不等式]] · [[27-热半群理论]] · [[24-弱极值原理与热方程初边值问题的唯一性]] · [[22-Cauchy问题与热核]]

[[24-弱极值原理与热方程初边值问题的唯一性]] 在有界区域上证明了唯一性，**没有用到任何增长假设** —— 弱极值原理把极值钉在抛物边界上，边界数据一给定，解就被锁死。

全空间不行。$\partial\Omega$ 消失了，柱面退到无穷远，剩下唯一能约束解的东西是初值；而热方程的热核 $H_t$ 对每个 $t>0$ 都有全支撑，"传播速度无穷"，初值根本管不住无穷远。于是"无穷远"变成一条**看不见但真实存在、上面没有任何数据**的边界。

后果不是技术性的：**唯一性真的会失效**。Tychonoff 1935 造出了 $\mathbb R^n_T$ 上齐次热方程的非零 $C^\infty$ 解，初值恒为 0。要救回唯一性，必须额外假设解在无穷远处的**增长**。

这一篇的路线：§26.1 说清"少了一条边"到底意味着什么，§26.2 造反例（Tychonoff），§26.3 解释反例从哪来（无穷远处的温差），§26.4 给出把反例排除掉的那条增长条件。

---

## 26.1 全空间少了一条边

**记号 26.1.1**　记

$$
\mathbb R^n_T:=\mathbb R^n\times(0,T],\qquad \mathbb R_T:=\mathbb R\times(0,T].
$$

有界区域的抛物边界（[[24-弱极值原理与热方程初边值问题的唯一性]] (24.1)）是**底盖加柱面**：

$$
\Gamma_T^*=\left(\overline\Omega\times\{0\}\right)\cup\left(\partial\Omega\times(0,T]\right).
$$

全空间时右边这一项没了，只剩底盖 $\mathbb R^n\times\{0\}$。弱极值原理的机制是"边界控制内部"，现在边界只剩下一条，控制力不够。

一个更本质的说法：**热方程使信息传播无穷快**。热核 $H_t$ 对任意小的 $t>0$ 都满足 $H_t>0$ 于全空间，所以 $t>0$ 时刻任何一点的值都依赖初值的**全部**数据，包括无穷远处。既然初值只给了"$t=0$ 时各点的值"，没给"无穷远处的边界值"，那信息就可以从无穷远自己灌进来 —— 合法的解不唯一。

这就是"适定性"里那三条（存在、唯一、连续依赖）中唯一性最容易失守的地方。

---

## 26.2 Tychonoff 反例：唯一性真的会失效

考察一维齐次问题

$$
\begin{cases}
\partial_tu=\partial_x^2u, & (x,t)\in\mathbb R_T,\\[2pt]
u(x,0)=0, & x\in\mathbb R .
\end{cases}\tag{26.1}
$$

平凡解 $u\equiv0$ 当然在。下面造一个非零的 $C^\infty$ 解。

### 26.2.1 一个在原点处「无限平坦」的函数

对 $\alpha>0$ 令

$$
u_0(t):=
\begin{cases}
e^{-1/t^{\alpha}}, & t>0,\\[2pt]
0, & t\le0 .
\end{cases}\tag{26.2}
$$

**引理 26.2.1**　$u_0\in C^\infty(\mathbb R)$，且 $u_0^{(k)}(0)=0$ 对一切 $k\ge0$。

**证明**　对 $t>0$ 归纳可写出

$$
u_0^{(k)}(t)=P_k\!\left(\frac1t\right)e^{-1/t^{\alpha}},
$$

其中 $P_k$ 是多项式（求导一次，指数上的 $-1/t^\alpha$ 贡献 $-1/t^\alpha$ 的导数即 $\alpha t^{-\alpha-1}$，与前面的多项式相乘仍是多项式乘指数）。令 $s=1/t\to+\infty$，则对任意整数 $m\ge0$ 有 $s^{m}e^{-s^{\alpha}}\to0$，故 $u_0^{(k)}(t)\to0$（$t\to0^+$）；在 $t<0$ 上 $u_0^{(k)}\equiv0$。于是 $u_0^{(k)}$ 在 $0$ 处连续且等于 $0$，对每个 $k$ 成立，故 $u_0\in C^\infty(\mathbb R)$ 且各阶导数在原点为零。$\square$

> **要点（"无限平坦"是整段构造的支点）**　$u_0$ 在 $0$ 处的 Taylor 级数是 $\sum_k u_0^{(k)}(0)t^k/k!=0$，但 $u_0\not\equiv0$ —— 即 **$u_0$ 在 $0$ 处不实解析**。正因如此，下面那个"以 $u_0$ 的各阶导数为系数"的级数才可能在 $t=0$ 处**假装成零**，而在 $t>0$ 处不为零。$C^\infty$ 与实解析之间的这条缝，就是反例的全部藏身之处。

### 26.2.2 级数解

令

$$
u(x,t):=\sum_{k=0}^{\infty}\frac{u_0^{(k)}(t)}{(2k)!}\,x^{2k}.\tag{26.3}
$$

**为什么是 $\dfrac{x^{2k}}{(2k)!}$**　热方程要求 $\partial_t$ 与 $\partial_x^2$ 把一项变成**同一族里的下一项**：

$$
\partial_x^2\left[\frac{x^{2k}}{(2k)!}\right]=\frac{x^{2k-2}}{(2k-2)!},\qquad \partial_t\left[u_0^{(k)}(t)\right]=u_0^{(k+1)}(t).
$$

两边各错开一位正好对上：空间求两次导把 $2k$ 降到 $2k-2$，时间求一次导把 $k$ 升到 $k+1$。

**验证方程。** 逐项求导（合法性见下），换指标 $j=k-1$：

$$
\partial_tu=\sum_{k\ge0}\frac{u_0^{(k+1)}(t)}{(2k)!}x^{2k},\qquad
\partial_x^2u=\sum_{k\ge1}\frac{u_0^{(k)}(t)}{(2k-2)!}x^{2k-2}=\sum_{j\ge0}\frac{u_0^{(j+1)}(t)}{(2j)!}x^{2j}.
$$

两式逐项相同，故 $\partial_tu=\partial_x^2u$ 于 $\mathbb R_T$。

**验证初值。** 由引理 26.2.1，$u(x,0)=\sum_k u_0^{(k)}(0)x^{2k}/(2k)!=0$。

**验证非零。** 取 $x=0$：

$$
u(0,t)=u_0(t)=e^{-1/t^{\alpha}}>0,\qquad t>0 .
$$

**验证收敛与逐项求导（讲稿整段省掉）。** 固定 $t>0$。$u_0$ 在圆盘 $|z-t|\le t/2$ 内解析且 $|u_0|\le1$，Cauchy 估计给出

$$
\left|u_0^{(k)}(t)\right|\le k!\left(\frac t2\right)^{-k}.
$$

又 $(2k)!\ge(k!)^2$，故 $\dfrac{k!}{(2k)!}\le\dfrac{1}{k!}$。于是

$$
\left|u(x,t)\right|\le\sum_{k\ge0}\left(\frac2t\right)^{k}\frac{k!\,x^{2k}}{(2k)!}\le\sum_{k\ge0}\frac{1}{k!}\left(\frac{2x^2}{t}\right)^{k}=\exp\left(\frac{2x^2}{t}\right).\tag{26.4}
$$

级数在 $\{t\ge t_0>0\}$ 的任意紧集上绝对一致收敛；同一估计对 $u_0^{(k+1)}$ 也成立，故逐项求导合法。$\square$

**结论**　(26.1) 有两个 $C^\infty$ 解：$u\equiv0$ 与 (26.3)。**唯一性失效。**

> **物理读法**　$t=0$ 时刻处处温度 $0$，$t>0$ 时刻在 $x=0$ 处凭空出现温度 $e^{-1/t^\alpha}$ 并向外扩散 —— 无源生热，反物理。讲稿的措辞是"$0$ 时刻恒温、$t>0$ 时刻热量却在流动"。这不是解的毛病，是"允许任意增长"这个前提的毛病。

### 26.2.3 它到底违反了哪一条

(26.4) 那个界是唯一的线索：指数常数是 $2/t$，随 $t\to0^+$ **无界**。所以 Tychonoff 解不在任何一个"$A$ 与 $t$ 无关"的增长类里：

$$
|u(x,t)|\le Ce^{A|x|^2}\qquad\text{对一切 }(x,t)\in\mathbb R^n_T\ \text{与同一组 }(A,C).\tag{26.5}
$$

（若取 $A=2/t_0$，只对 $t\ge t_0$ 成立；反例恰恰在 $t\to0^+$ 处越界。）

对比真正由卷积给出的解：$\varphi$ 有紧支集时 $|u(x,t)|\sim e^{-C|x|^2}$（高斯**衰减**，见订正 2），$\varphi$ 仅有界时 $|u(x,t)|\le\lVert\varphi\rVert_{L^\infty}$（有界）。两者都远小于 $e^{2x^2/t}$。**多出来的那个解藏在无穷远处**，它的增长恰好落在"所有固定指数的界都管不到、但又不是完全无界"的那条缝里。

§26.4 的定理说的正是：**把 $A$ 钉死，唯一性就回来。**

---

## 26.3 为什么会有反例：无穷远处的温差

反例不是凭空的。既然 $t>0$ 时每一点都感受得到无穷远，那就该问：**无穷远处的温度到底给了多少？** 把它量化，就得到下面这套障碍函数论证。

### 26.3.1 温差 $b_R$

设 $w$ 解

$$
\begin{cases}
\partial_tw=\Delta w, & (x,t)\in\mathbb R^n_T,\\[2pt]
w(x,0)=0, & x\in\mathbb R^n .
\end{cases}
$$

对 $R>0$ 令 $b_R$ 为下列**有界区域**问题的解：

$$
\begin{cases}
\partial_tb_R=\Delta b_R, & (x,t)\in B_R\times(0,T],\\
b_R(x,0)=0, & x\in B_R,\\
b_R(x,t)=\lvert w(x,t)\rvert, & x\in\partial B_R,\ t\in[0,T].
\end{cases}\tag{26.6}
$$

$b_R$ 称为 $\partial B_R$ 上的**温差**：它只被球面上那一圈数据驱动，度量的是"球面温度对球内的最大可能影响"。把 $\lim_{R\to\infty}b_R$ 理解为 Cauchy 问题**无穷远处的温差**。

**引理 26.3.1**　对任意 $R>0$，

$$
0\le\lvert w(x,t)\rvert\le b_R(x,t)\qquad\text{于 }B_R\times[0,T].
$$

**证明**　

$b_R\ge0$：初值与边值都非负，弱极值原理（[[24-弱极值原理与热方程初边值问题的唯一性]] 定理 24.2.1）。

$\lvert w\rvert\le b_R$：$w$ 与 $-w$ 都解热方程，从而是（弱）**下解**；**两个下解取 $\max$ 仍是下解**，故 $\lvert w\rvert=\max\{w,-w\}$ 是下解。$b_R$ 是解、从而是上解，且在抛物边界上 $b_R\ge\lvert w\rvert$（侧边取等，底盖两边都是 $0$）。比较原理给出 $\lvert w\rvert\le b_R$。$\square$

> **要点（这一步不能跳）**　讲稿直接从"$0\le b_R\le\Phi_R$"跳到"$w\equiv0$"，中间缺的正是 $\lvert w\rvert\le b_R$。没有它，$b_R\to0$ 对 $w$ 什么也说明不了。而它成立的理由不是"$w$ 是解"（$\lvert w\rvert$ 一般**不解**热方程，在 $w$ 变号处会多出 $-\lvert\nabla w\rvert$ 型的负源），而是"$\max$ 保持下解性"这条结构性事实 —— 这是极值原理类论证里最常用、也最容易被跳过的一步。

### 26.3.2 障碍函数 $\Phi_R$ 是怎么想出来的

要让 $b_R\to0$，只需给 $b_R$ 找一个随 $R\to\infty$ 消失的上界。

**第一步：找"增长类对应的那个解"。** 设 $a>0$，记 $\tau:=1-4at$，令

$$
\psi_a(x,t):=\tau^{-p}\exp\left(\frac{a\lvert x\rvert^2}{\tau}\right)\qquad(\tau>0).
$$

直接算（$r^2=\lvert x\rvert^2$，注意 $\partial_t\tau=-4a$）：

$$
\partial_t\psi_a=\psi_a\left[4ap\,\tau^{-1}+4a^2r^2\tau^{-2}\right],\qquad
\Delta\psi_a=\psi_a\left[2an\,\tau^{-1}+4a^2r^2\tau^{-2}\right].
$$

两式相减，$r^2$ 项完全对消：

$$
\left(\partial_t-\Delta\right)\psi_a=\frac{2a\,(2p-n)}{\tau}\,\psi_a .\tag{26.7}
$$

**所以 $\psi_a$ 解热方程当且仅当 $p=n/2$。**（$n=1$ 时正是 $p=1/2$ —— 讲稿写的那个数只在 $n=1$ 对，见订正 3。）

**$\psi_a$ 是什么？** 看 $t=0$：$\psi_a(x,0)=e^{a\lvert x\rvert^2}$ —— **它就是增长类边界函数本身**。$\psi_a$ 是它的向前解析延拓，在 $t\nearrow1/(4a)$ 时爆破。所以"$e^{a\lvert x\rvert^2}$ 型的初值只能存活到 $t=1/(4a)$"这件事，在 $\psi_a$ 上一个字就看清了：**增长条件越宽（$a$ 越大），能活的时间越短。**

**第二步：造障碍。** 取 $a=A+\varepsilon$（$\varepsilon>0$ 是留给收缩因子的余量），令

$$
\Phi_R(x,t):=Ce^{-\varepsilon R^2}\psi_{A+\varepsilon}(x,t)
=Ce^{-\varepsilon R^2}\left(1-4(A+\varepsilon)t\right)^{-n/2}\exp\left(\frac{(A+\varepsilon)\lvert x\rvert^2}{1-4(A+\varepsilon)t}\right).\tag{26.8}
$$

前提是 $0\le t\le T<\dfrac{1}{4(A+\varepsilon)}$，即 $\tau\in(0,1]$（讲稿未写，见订正 4）。在此前提下三条成立：

1. $\partial_t\Phi_R=\Delta\Phi_R$ —— 由 (26.7) 与 $p=n/2$；
2. 底盖：$\Phi_R(x,0)=Ce^{-\varepsilon R^2}e^{(A+\varepsilon)\lvert x\rvert^2}\ge0=b_R(x,0)$；
3. 柱面：$\lvert x\rvert=R$ 且 $\tau\le1$ 给出 $\dfrac{(A+\varepsilon)R^2}{\tau}\ge(A+\varepsilon)R^2$，故

$$
\Phi_R\big|_{\partial B_R\times[0,T]}\ \ge\ Ce^{-\varepsilon R^2}e^{(A+\varepsilon)R^2}=Ce^{AR^2}\ \ge\ \lvert w\rvert .
$$

比较原理 ⟹

$$
0\le b_R\le\Phi_R\qquad\text{于 }B_R\times[0,T].\tag{26.9}
$$

**第三步：令 $R\to\infty$。** 固定 $(x,t)$，$\Phi_R(x,t)=Ce^{-\varepsilon R^2}\cdot[\text{与 }R\text{ 无关的因子}]\to0$。于是

$$
w(x,t)\le b_R(x,t)\le\Phi_R(x,t)\to0\quad(R\to\infty),
$$

对 $-w$ 同法得 $w(x,t)\ge0$。故 $w\equiv0$。

> **要点（"$\lim_R\Phi_R\equiv0$"与"$\Phi_R|_{\partial B_R}\ge Ce^{AR^2}$"为什么不矛盾）**　两处用的是**完全不同的 $R$ 依赖**：
>
> - 球面上 $\lvert x\rvert=R$ 跟着 $R$ 一起跑，指数是 $\frac{(A+\varepsilon)R^2}{\tau}-\varepsilon R^2\ge AR^2$，**增长**；
> - 固定内部点上 $\lvert x\rvert$ 不动，指数只剩 $-\varepsilon R^2$，**衰减**。
>
> 所以 $e^{-\varepsilon R^2}$ 不是"多余的小因子"，而是整个论证的开关：拿掉它，球面上的下界还在，但内部的界就永远是 $e^{a\lvert x\rvert^2}$ 而不是 $0$。$A$ 与 $\varepsilon$ 各管一头 —— **$A$ 定球面高度（必须压过 $\lvert w\rvert$），$\varepsilon$ 定收缩速度（必须让内部趋于 $0$）**。这也是唯一可以用到"$A$ 严格固定"的地方：若 $A$ 依赖 $t$，$\Phi_R$ 就写不出来。

---

## 26.4 具有增长控制的最大值原理

**定理 26.4.1**　设 $u$ 是 $\mathbb R^n_T$ 上热方程的解，且存在 $A>0$、$C>0$ 使

$$
u(x,t)\le Ce^{A\lvert x\rvert^2},\qquad\forall (x,t)\in\mathbb R^n_T .
$$

则

$$
\sup_{\mathbb R^n_T}u\ \le\ \sup_{x\in\mathbb R^n}u(x,0).\tag{26.10}
$$

注意只有**上**界。这正是"上解"框架需要的形式：给定零初值时推出 $u\le0$ 于全时空。双侧的估计要分别对 $u$ 与 $-u$ 各用一次（§26.4.4）。

### 26.4.1 障碍函数的选法：把热核倒着走

取

$$
\sigma:=2T-t\in[T,\,2T],\qquad
g(x,t):=\frac{\delta}{\sigma^{n/2}}\exp\left(\frac{\lvert x\rvert^2}{4\sigma}\right),\qquad
v_\delta:=u-g\le u .\tag{26.11}
$$

**为什么 $g$ 解热方程。** $g$ 就是**热核在"负时间" $s=t-2T<0$ 处的取值**：把

$$
H(x,s)=(4\pi s)^{-n/2}e^{-\frac{\lvert x\rvert^2}{4s}}
$$

里的 $s$ 换成 $t-2T=-\sigma$，指数上的负号与 $s<0$ 相消，就得到 $+\frac{\lvert x\rvert^2}{4\sigma}$。$H(\cdot,s)$ 对 $s\ne0$ 都解热方程（公式本身对 $s<0$ 也成立），而 $t\mapsto t-2T$ 的导数是 $1$，所以 $g$ 对 $t$ 也解热方程。

**为什么用 $2T$ 而不是 $T$。** $s=0$ 是唯一的奇点。取 $s=t-2T$ 时它落在 $t=2T$，在 $[0,T]$ 之外，全段安全；若取 $s=t-T$，奇点正好落在右端点 $t=T$，$(0,T]$ 上就无定义。

**逐项验算（讲稿只列结果）。** 由 $\partial_t\sigma=-1$：

$$
\partial_t\left(\sigma^{-n/2}\right)=\frac n2\sigma^{-n/2-1},\qquad
\partial_t\left(\frac{\lvert x\rvert^2}{4\sigma}\right)=\frac{\lvert x\rvert^2}{4\sigma^2},
$$

$$
\partial_tg=\delta e^{\frac{\lvert x\rvert^2}{4\sigma}}\left[\frac n2\sigma^{-n/2-1}+\frac{\lvert x\rvert^2}{4\sigma^2}\sigma^{-n/2}\right];
$$

$$
\partial_{x_i}g=\delta\sigma^{-n/2}e^{\frac{\lvert x\rvert^2}{4\sigma}}\frac{x_i}{2\sigma},\qquad
\partial^2_{x_ix_i}g=\delta\sigma^{-n/2}e^{\frac{\lvert x\rvert^2}{4\sigma}}\left[\frac{1}{2\sigma}+\frac{x_i^2}{4\sigma^2}\right],
$$

对 $i$ 求和（$\sum_i x_i^2=\lvert x\rvert^2$）：

$$
\Delta g=\delta e^{\frac{\lvert x\rvert^2}{4\sigma}}\left[\frac n2\sigma^{-n/2-1}+\frac{\lvert x\rvert^2}{4\sigma^2}\sigma^{-n/2}\right]=\partial_tg .
$$

两处符号的相消值得单独记：$\sigma$ 的**时间反演**（$\partial_t\sigma=-1$）与指数上的**正号**（$+\frac{\lvert x\rvert^2}{4\sigma}$）各出一个负号，乘起来抵消。**于是 $v_\delta$ 也解热方程**：$(\partial_t-\Delta)v_\delta=0$。

### 26.4.2 抛物边界上的两个估计

在 $B_R\times(0,T]$ 上用弱极值原理：

$$
\sup_{B_R\times[0,T]}v_\delta\ \le\ \sup_{\Gamma_T}v_\delta ,\qquad
\Gamma_T=\left(\overline{B_R}\times\{0\}\right)\cup\left(\partial B_R\times[0,T]\right).
$$

**底盖（$t=0$，$\lvert x\rvert\le R$）。** 此时 $\sigma=2T$，

$$
v_\delta(x,0)=u(x,0)-\frac{\delta}{(2T)^{n/2}}e^{\frac{\lvert x\rvert^2}{8T}}\le u(x,0)\le\sup_{x\in\mathbb R^n}u(x,0)=:M .
$$

**柱面（$\lvert x\rvert=R$，$0<t\le T$）。** 由 $\sigma\le2T$ 同时得 $\sigma^{-n/2}\ge(2T)^{-n/2}$ 与 $\frac{R^2}{4\sigma}\ge\frac{R^2}{8T}$，故

$$
v_\delta\ \le\ Ce^{AR^2}-\frac{\delta}{(2T)^{n/2}}e^{\frac{R^2}{8T}}
=e^{AR^2}\left[C-\frac{\delta}{(2T)^{n/2}}e^{R^2\left(\frac{1}{8T}-A\right)}\right].
$$

若 $\frac{1}{8T}>A$，即 $T<\frac{1}{8A}$，方括号随 $R\to+\infty$ 趋于 $-\infty$，故对充分大的 $R$ 有 $v_\delta\le0\le M$ 于柱面。

于是对该 $R$ 有 $\sup_{B_R\times[0,T]}v_\delta\le M$。固定 $(x,t)$，取足够大的 $R>\lvert x\rvert$ 得 $v_\delta(x,t)\le M$；再由 $v_\delta(x,t)\to u(x,t)$（$\delta\to0^+$，**逐点**收敛）得 $u(x,t)\le M$。$\square$（此段只证了 $T<\frac{1}{8A}$ 的情形）

### 26.4.3 去掉 $T<1/(8A)$

$T<\frac{1}{8A}$ 不是结论的一部分，只是**障碍函数的寿命限制** —— §26.3.2 已经说明增长类 $e^{a\lvert x\rvert^2}$ 只能活到 $t=1/(4a)$，这里 $\frac{1}{8A}<\frac{1}{4A}$，是同一个现象。

取 $T_0:=\frac{1}{16A}<\frac{1}{8A}$，正整数 $N$ 使 $(N-1)T_0<T\le NT_0$，把 $[0,T]$ 切成

$$
[0,T_0]\cup[T_0,2T_0]\cup\cdots\cup[(N-1)T_0,\,T].
$$

每段长度 $\le T_0<\frac{1}{8A}$。在第 $i$ 段上作时间平移

$$
\tilde u(x,\tau):=u(x,\,iT_0+\tau),\qquad \tau\in(0,T_i],\quad T_i:=\min\{(i+1)T_0,T\}-iT_0\le T_0 .
$$

$\tilde u$ 仍解热方程，且仍满足**同一个**增长界 $\tilde u\le Ce^{A\lvert x\rvert^2}$，故上面的结论对它成立：

$$
\sup_{x\in\mathbb R^n,\ iT_0\le t\le\min\{(i+1)T_0,T\}}u(x,t)\ \le\ \sup_{x\in\mathbb R^n}u(x,\,iT_0).
$$

逐段串起来即得 (26.10)。$\square$

> **要点（分段法能成立的全部理由）**　**增长常数 $A$ 与 $t$ 无关。** 每一段平移后面对的都是同一个 $A$，所以每段都能独立地取 $T_i\le T_0<\frac{1}{8A}$。若增长界写成 $e^{A(t)\lvert x\rvert^2}$ 而 $A(t)$ 随时间变（Tychonoff 反例正是 $A(t)=2/t$，$t\to0^+$ 时爆掉），分段法立刻失效 —— 而反例也正因此存在。**"$A$ 固定"这一条，就是定理与反例之间的那条线。**

> **要点（为什么是 $\frac{1}{16A}$ 而不是 $\frac{1}{8A}$）**　$T_0$ 必须**严格**小于 $\frac{1}{8A}$，因为最后一段 $[(N-1)T_0,T]$ 的长度可能正好等于 $T_0$。取一半是最省事的写法。取 $\frac{1}{16A}$ 而非 $\frac{1}{10A}$ 纯属整齐。

### 26.4.4 推论：增长类中的唯一性

**推论 26.4.2**　设 $f\in C(\mathbb R^n_T)$、$g\in C(\mathbb R^n)$，则 Cauchy 问题

$$
\begin{cases}
\partial_tu-\Delta u=f, & (x,t)\in\mathbb R^n_T,\\[2pt]
u(x,0)=g(x), & x\in\mathbb R^n
\end{cases}\tag{26.12}
$$

在增长类 $\left\{u:\ \lvert u(x,t)\rvert\le Ce^{A\lvert x\rvert^2}\right\}$（$C,A$ 固定）中**至多有一个**解。

**证明**　设 $u_1,u_2$ 是两个解，$w:=u_1-u_2$。则 $w$ 解齐次方程、$w(\cdot,0)=0$，且 $\lvert w\rvert\le2Ce^{A\lvert x\rvert^2}$。对 $w$ 用定理 26.4.1：

$$
\sup_{\mathbb R^n_T}w\le\sup_{x\in\mathbb R^n}w(x,0)=0 ;
$$

对 $-w$ 用同一定理（它满足同样的增长界）：

$$
\sup_{\mathbb R^n_T}(-w)\le0\quad\Longleftrightarrow\quad w\ge0 .
$$

故 $w\equiv0$。$\square$

> **要点**　定理 26.4.1 是**单侧**的，唯一性必须靠"两次单侧"凑出双侧 —— 这是所有极值原理型唯一性证明的共同套路（对照 [[24-弱极值原理与热方程初边值问题的唯一性]] 推论 2）。也正因如此，定理的假设只写了 $u\le Ce^{A\lvert x\rvert^2}$ 而没有 $\lvert u\rvert$：多写一个绝对值，就等于把两次论证偷偷塞进一次里，反而看不出哪一步用了单侧性。

---

## 26.5 三类增长条件

| 增长假设 | 唯一性 | 依据 |
|---|---|---|
| $\lvert u\rvert\le M$（有界） | 成立 | 有界 ⟹ $\lvert u\rvert\le Me^{A\lvert x\rvert^2}$ 对任意 $A$ |
| $\lvert u\rvert\le Ce^{A\lvert x\rvert^2}$，$A,C$ 与 $t$ 无关 | 成立 | 定理 26.4.1 + 推论 26.4.2 |
| $\lvert u\rvert\le Ce^{a(t)\lvert x\rvert^2}$，$a(t)\to+\infty$（$t\to0^+$） | **失效** | Tychonoff 反例，$a(t)=2/t$ |
| 无任何增长假设 | **失效** | 同上 |

一句话：**全空间热方程 Cauchy 问题在 $C^\infty$ 类里不适定，只在"指数二次型、常数钉死"的增长类里适定。**

对照 Laplace 方程：椭圆问题没有时间方向，也就没有"从无穷远单向灌进来的信息**在未给定边界数据时**"，所以 Tychonoff 型现象在椭圆情形不存在（[[13-奇点可去性与局部解析性]] 的可去奇点定理是另一回事，它讲的是"点奇点"，不是"无穷远"）。

---

## 26.6 与有界区域的对照

| | 有界区域（[[24-弱极值原理与热方程初边值问题的唯一性]]） | 全空间（本篇） |
|---|---|---|
| 抛物边界 | 底盖 $\cup$ 柱面 | 只剩底盖；"无穷远"是隐形边界，无数据 |
| 弱极值原理 | 无条件成立 | 需要增长控制，否则反例 |
| 唯一性 | 无条件 | 增长类内唯一 |
| 控制机制 | 侧边数据把极值钉在边界 | 无穷远温差被障碍函数压成 $0$ |
| 障碍函数 | $\varepsilon x_1^2$（多项式，永远正） | $\Phi_R$（指数二次型，**寿命有限**） |
| 反例 | 无 | Tychonoff |

最后一行是全部差别的根源：**有界区域用的障碍函数 $\varepsilon x_1^2$ 对任意大的 $t$ 都存在；全空间的障碍函数 $\Phi_R$ 活不过 $t=1/(4a)$。** 因为全空间没有侧边可以"截断"，只能靠函数自己在无穷远处压低 —— 而"在无穷远压低"这件事，指数二次型只能维持有限时间。

---

## 订正

1. **§4.1 的级数解只写了"则 $u_t=u_{xx}$、$u(x,0)=0$"，收敛性与逐项求导的合法性完全没验。** 这不是可以省的技术细节：级数 (26.3) 的系数是 $u_0$ 的 $k$ 阶导数，若不先给 $u_0^{(k)}(t)$ 的增长界，$\sum u_0^{(k)}(t)x^{2k}/(2k)!$ 连逐点收敛都不显然。本篇 §26.2.2 用 Cauchy 估计补完，并顺手得到 (26.4) —— 而这个界恰好就是 §26.2.3 判定"反例违反哪条增长条件"的唯一依据。**验收敛不是补手续，是反例论证的另一半。**

2. **§4.1 末句"若 $\varphi$ 有界、连续，则对固定的 $t>0$，$\lvert u(x,t)\rvert\sim e^{-C\lvert x\rvert^2}$"不成立。** 取 $\varphi\equiv1$ 则 $u\equiv1$，毫无衰减。正确说法：$\varphi$ 有**紧支集**（或有足够衰减）时才有 $t$ 固定、$\lvert x\rvert\to\infty$ 时的高斯衰减；$\varphi$ 仅有界连续时只能得 $\lvert u\rvert\le\lVert\varphi\rVert_{L^\infty}$。好在这句话真正要用的只是"真解落在增长类 (26.5) 里"，而这一条对 $\varphi$ 仅有界就成立（有界自然满足指数二次型的界）。

3. **§4.2 的 $\Phi_R$ 把指数写成 $1/2$，在 $n$ 维下是错的，应为 $n/2$。** 由 (26.7)，$\psi_a$ 解热方程的**充要条件**是 $p=n/2$。取 $p=1/2$、$n\ge2$ 时 $(\partial_t-\Delta)\Phi_R=\frac{2a(1-n)}{\tau}\Phi_R<0$ —— 是**下解**，比较原理的方向反了，$b_R\le\Phi_R$ 推不出来。所以这不是"差个常数、无伤大雅"：$1/2$ 只在 $n=1$ 恰好对，而这个例子是在 $\mathbb R^n$ 里用的。判据：**先算 $(\partial_t-\Delta)\psi_a$，看幂次能不能把 $\partial_t$ 贡献的 $4ap$ 与 $\Delta$ 贡献的 $2an$ 配平。**

4. **§4.2 缺一个前提：$T\le\dfrac{1}{4(A+\varepsilon)}$。** 没有它，$\tau=1-4(A+\varepsilon)t$ 可能变号，$\Phi_R$ 在 $\tau<0$ 时指数上的 $\frac{a\lvert x\rvert^2}{\tau}$ 变负、柱面上的下界 $\Phi_R\ge Ce^{AR^2}$ 直接崩掉。这个前提与 §4.3 的 $T<\frac{1}{8A}$ 是同一现象的两处表现 —— 增长类的"陪伴解"只能在 $0\le t<1/(4a)$ 上存在（§26.3.2）。既然 §4 最后要用分段把它去掉，§4.2 这里至少该点明"先设 $T$ 小"。

5. **"手稿中记 $\Phi_R$ 满足热方程"——"记"字说明这是未核实的转述。** 这句话是全节的支点（不成立则比较原理无从下手），却只以"手稿"作保。本篇 §26.3.2 给了完整计算，并把 $p=n/2$ 这个充要条件单独写出 —— 顺带解释了"为什么恰好是这个形式"（它是增长类边界函数 $e^{a\lvert x\rvert^2}$ 的向前解析延拓）。

6. **从 $0\le b_R\le\Phi_R$ 到 $w\equiv0$ 中间缺了 $\lvert w\rvert\le b_R$。** 见引理 26.3.1 及其要点。缺这一步时，"$b_R\to0$"与"$w\equiv0$"之间没有逻辑通道。另外 $b_R\ge0$ 也需要一句（初值边值非负 + 弱极值原理）。

7. **§4.3 的 $v_\delta$ 定义后直接写 $v_\delta\le u$，没说明为什么。** 理由是 $g>0$（$\sigma>0$、指数为正、$\delta>0$）。看似显然，但这个不等号方向是后面"$v_\delta\le M$ ⟹ $u\le M$"唯一的依据，写清楚只需半句。

8. **§4.3 的分段只说"平移时间变量、按该段长度定义辅助函数"，没有点明分段法成立的命门是"$A$ 与 $t$ 无关"。** 如果读者以为分段是纯技巧，就会误以为任何增长界都能这样处理 —— 而 Tychonoff 反例恰恰是 $A(t)=2/t$、分段失效。**"$A$ 固定"才是定理与反例的分界线**，这句话必须在分段处说出来（本篇 §26.4.3 要点）。

9. **$\frac{1}{8A}$ 与 $T_0=\frac{1}{16A}$ 两个常数的关系没解释。** 为什么不是 $T_0=\frac{1}{8A}$？因为最后一段 $[(N-1)T_0,T]$ 可能正好取满 $T_0$，需要**严格**小于。取一半是最省事的写法。

---

系列：[[25-热方程的梯度估计与Harnack不等式]] · [[27-热半群理论]] · [[24-弱极值原理与热方程初边值问题的唯一性]] · [[22-Cauchy问题与热核]]
