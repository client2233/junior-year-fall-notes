# 25 热方程的梯度估计与 Harnack 不等式

> 系列：[[24-弱极值原理与热方程初边值问题的唯一性]] · [[26-Cauchy问题的唯一性与Tychonoff反例]] · [[04-调和函数·Harnack不等式与内部梯度估计]] · [[22-Cauchy问题与热核]]

这一篇是 [[04-调和函数·Harnack不等式与内部梯度估计]] 的抛物版。04 篇对调和函数做了两件事：**内部梯度估计**（$\lvert\nabla\log u\rvert\le C/R$）与 **Harnack 不等式**（$u(x)\le Cu(y)$）。这里对热方程做同样两件事，结论形状几乎平行，但**每一项都多一个时间导数**，常数的位置也从 $R$ 换成了 $t$：

$$
\lvert\nabla\log u\rvert^2-\partial_t\log u\le\frac{n}{2t}.\tag{25.1}
$$

两个证明路线值得对照着看：

- **§25.1 用热核表示 + 加权 Cauchy–Schwarz** —— 只要初值非负，一行 Cauchy–Schwarz 就完事，且常数 $n/2t$ 显然最优；
- **§25.2 用 Bochner 恒等式 + 极值原理** —— 不需要初值，只要解**正**，但计算量大得多。

两者给出**同一个不等式**。这种「表示法给结论、Bochner 法给一般性」的分工，在 Laplace 情形就是 [[04-调和函数·Harnack不等式与内部梯度估计]] 里 Poisson 核证法与 Bochner 证法的关系。

---

## 25.1 全局梯度估计（初值非负）

**定理 25.1.1**　设 $\varphi\in C_b(\mathbb R^n)$、$\varphi\ge0$、$\varphi\not\equiv0$，令

$$
u(x,t)=\left(H_t*\varphi\right)(x)=\frac{1}{(4\pi t)^{n/2}}\int_{\mathbb R^n}e^{-\frac{|x-y|^2}{4t}}\varphi(y)\,dy,\qquad t>0,\tag{25.2}
$$

其中 $H_t(x):=H(x,t)=(4\pi t)^{-n/2}e^{-\frac{|x|^2}{4t}}$（[[22-Cauchy问题与热核]] 定理 (22.7)）。则 $u(x,t)>0$，且对一切 $(x,t)\in\mathbb R^n\times(0,+\infty)$，

$$
\boxed{\ \frac{\left|\nabla u\right|^2}{u^2}-\frac{\partial_tu}{u}\le\frac{n}{2t}\ }\tag{25.3}
$$

### 25.1.1 先证 $u>0$

取 $y_0$ 使 $\varphi(y_0)>0$，由 $\varphi$ 连续存在 $r>0$ 使 $\varphi\ge\varphi(y_0)/2$ 于 $B_r(y_0)$。被积函数非负，故只需在 $B_r(y_0)$ 上估下界：

$$
u(x,t)\ge\frac{1}{(4\pi t)^{n/2}}\int_{B_r(y_0)}e^{-\frac{|x-y|^2}{4t}}\varphi(y)\,dy\ge\frac{\varphi(y_0)}{2}\cdot\frac{\left|B_r\right|}{(4\pi t)^{n/2}}e^{-\frac{\left(|x-y_0|+r\right)^2}{4t}}>0 .
$$

（指数上的 $(|x-y_0|+r)^2$ 来自三角不等式 $|x-y|\le|x-y_0|+r$。）**这一步是必需的**：(25.3) 两边都除以 $u$，$u$ 必须有下界。

### 25.1.2 热核自身的恒等式

**引理 25.1.2**　对 $t>0$、$z\in\mathbb R^n$，

$$
\frac{\left|\nabla H(z,t)\right|^2}{H(z,t)^2}-\frac{\partial_tH(z,t)}{H(z,t)}=\frac{n}{2t}.\tag{25.4}
$$

**证明**　由 $H=(4\pi t)^{-n/2}e^{-\frac{|z|^2}{4t}}$ 直接算（与 [[22-Cauchy问题与热核]] §22.5 第 3 条同一套计算）：

$$
\nabla H=H\cdot\left(-\frac{z}{2t}\right)\ \Longrightarrow\ \frac{\left|\nabla H\right|^2}{H^2}=\frac{|z|^2}{4t^2},
$$

$$
\partial_tH=H\left(-\frac{n}{2t}+\frac{|z|^2}{4t^2}\right)\ \Longrightarrow\ \frac{\partial_tH}{H}=-\frac{n}{2t}+\frac{|z|^2}{4t^2}.
$$

两式相减，$|z|^2/4t^2$ 对消，得 $n/2t$。$\square$

> **要点**　(25.4) 说明**热核本身取到 (25.3) 的等号**。也就是说 $n/2t$ 这个常数不能再改进 —— 只要初值"足够集中"（$\varphi\to\delta_{y_0}$），$u$ 就趋近于一个平移后的热核，(25.3) 的比值趋于 $n/2t$。**估计的最优性不用另证，看引理 25.1.2 就够了。**

### 25.1.3 加权 Cauchy–Schwarz

把 (25.3) 左端写成 $\log u$ 的形式会更清楚：

$$
\frac{\left|\nabla u\right|^2}{u^2}-\frac{\partial_tu}{u}=\left|\nabla\log u\right|^2-\partial_t\log u .
$$

关键是把 $\dfrac{\nabla u}{u}$ 写成**一个概率测度下的平均**。由 (25.2)，

$$
\nabla u(x,t)=\int_{\mathbb R^n}\nabla H(x-y,t)\varphi(y)\,dy,\qquad u(x,t)=\int_{\mathbb R^n}H(x-y,t)\varphi(y)\,dy .
$$

两式相除，并在被积函数里乘除同一个 $H$：

$$
\frac{\nabla u}{u}=\int_{\mathbb R^n}\frac{\nabla H(x-y,t)}{H(x-y,t)}\cdot\underbrace{\frac{H(x-y,t)\varphi(y)}{u(x,t)}}_{=:\ \nu_{x,y,t}}dy .\tag{25.5}
$$

**$\nu_{x,y,t}$ 是概率密度**：由 $\varphi\ge0$、$H>0$、$u>0$ 得 $\nu\ge0$，且

$$
\int_{\mathbb R^n}\nu_{x,y,t}\,dy=\frac{1}{u(x,t)}\int_{\mathbb R^n}H(x-y,t)\varphi(y)\,dy=\frac{u(x,t)}{u(x,t)}=1 .\tag{25.6}
$$

于是对 (25.5) 用 Cauchy–Schwarz（等价于 Jensen 不等式：$\left(\int g\,d\nu\right)^2\le\int g^2\,d\nu$），取 $g=\nabla H/H$：

$$
\frac{\left|\nabla u\right|^2}{u^2}\le\int_{\mathbb R^n}\frac{\left|\nabla H(x-y,t)\right|^2}{H(x-y,t)^2}\,\nu_{x,y,t}\,dy .
$$

代 (25.4)（把 $z$ 换成 $x-y$）：

$$
\frac{\left|\nabla u\right|^2}{u^2}\le\int_{\mathbb R^n}\left(\frac{\partial_tH(x-y,t)}{H(x-y,t)}+\frac{n}{2t}\right)\nu_{x,y,t}\,dy=\frac{1}{u}\int_{\mathbb R^n}\partial_tH(x-y,t)\varphi(y)\,dy+\frac{n}{2t}\cdot\underbrace{\int\nu\,dy}_{=1},
$$

末步用了 $\nu$ 的定义把 $H$ 约掉、以及 (25.6)。而 $\int\partial_tH\,\varphi=\partial_tu$（$t>0$ 时导数可搬进积分号，理由同 [[22-Cauchy问题与热核]] 命题 22.6.1），故

$$
\frac{\left|\nabla u\right|^2}{u^2}\le\frac{\partial_tu}{u}+\frac{n}{2t},
$$

移项即 (25.3)。$\square$

> **要点（这个证明的全部内容）**　只有三步：**把梯度写成概率平均（(25.5)+(25.6)）→ Cauchy–Schwarz → 代热核自身的等式 (25.4)**。没有任何"技巧"，全靠 $\nu$ 的归一化。**记住这个套路**：凡是能写成"某个核的加权平均"的量，梯度估计就是一次 Cauchy–Schwarz。
>
> 对照 [[04-调和函数·Harnack不等式与内部梯度估计]]：那里对 Poisson 核做同样的事，得到 $\lvert\nabla u\rvert/u\le C/R$；这里对热核做，得到的是**平方**形式的 (25.3)。差别在于 Poisson 核的 $\lvert\nabla H\rvert/H$ 只被 $C/R$ 控制（阶数恰好），热核的 $\lvert\nabla H\rvert/H=|z|/2t$ 无界，所以只能控制平方平均 —— 而平方平均恰好够用。

---

## 25.2 一般正解：Bochner 论证

讲稿在 §25.1 之后加了一句：「若 $u(x,t)$ 为 $\mathbb R^n\times(0,T]$ 上**一般的正解**，不考虑初值，同样有梯度估计。」然后直接陈述结论，**没有证明**。这里补上 —— 因为它是这一篇最漂亮的一段计算，而且它说明 (25.3) **不需要初值非负**，只需要**解正**。

**定理 25.2.1**　设 $u>0$ 于 $\mathbb R^n\times(0,T]$ 且 $\partial_tu=\Delta u$。则对一切 $(x,t)\in\mathbb R^n\times(0,T]$，

$$
\left|\nabla\log u\right|^2-\partial_t\log u\le\frac{n}{2t}.\tag{25.7}
$$

**证明**　记 $f:=\log u$。

**第一步：$f$ 满足的方程。** 由 $u=e^f$、$\Delta u=e^f\left(\Delta f+\lvert\nabla f\rvert^2\right)$，

$$
f_t=\frac{u_t}{u}=\frac{\Delta u}{u}=\Delta f+\left|\nabla f\right|^2\quad\Longleftrightarrow\quad \left|\nabla f\right|^2-f_t=-\Delta f .\tag{25.8}
$$

**所以 (25.7) 等价于 $-\Delta f\le\dfrac{n}{2t}$ 于 $\mathbb R^n\times(0,T]$。** 记住这一点：整段论证的目标就是给 $-\Delta\log u$ 一个上界。

**第二步：构造待估量并算它的演化。** 令

$$
F(x,t):=t\left(\left|\nabla f\right|^2-f_t\right)=-t\,\Delta f .
$$

用 $\Delta\lvert\nabla f\rvert^2=2\lvert\nabla^2f\rvert^2+2\langle\nabla f,\nabla\Delta f\rangle$（Bochner 恒等式）与 $(\partial_t-\Delta)f=\lvert\nabla f\rvert^2$，直接计算给出

$$
\left(\partial_t-\Delta\right)F=\frac{F}{t}-2t\left|\nabla^2f\right|^2+2\left\langle\nabla f,\nabla F\right\rangle .\tag{25.9}
$$

**第三步：用 $\lvert\nabla^2f\rvert^2\ge\frac1n(\Delta f)^2$。** 由 $\Delta f=-F/t$，

$$
\left|\nabla^2f\right|^2\ge\frac{1}{n}\left(\Delta f\right)^2=\frac{F^2}{nt^2}\quad\Longrightarrow\quad-2t\left|\nabla^2f\right|^2\le-\frac{2F^2}{nt}.
$$

代入 (25.9) 并移项：

$$
\left(\partial_t-\Delta-2\left\langle\nabla f,\nabla\cdot\right\rangle\right)F\ \le\ \frac{F}{t}-\frac{2F^2}{nt}.\tag{25.10}
$$

**第四步：极值原理。** 设 $F$ 在内部点取到**正的**最大值。则该点处 $\nabla F=0$、$\Delta F\le0$、$\partial_tF\ge0$，于是 (25.10) 左端 $\ge0$，故

$$
0\le\frac{F}{t}-\frac{2F^2}{nt}\quad\Longrightarrow\quad F\le\frac{n}{2}.
$$

即 $t\left(\lvert\nabla f\rvert^2-f_t\right)\le\frac n2$，除以 $t>0$ 得 (25.7)。$\square$

> **要点（(25.10) 里那个一阶项是怎么来的）**　注意左端不是 $\partial_t-\Delta$ 而是 $\partial_t-\Delta-2\langle\nabla f,\nabla\cdot\rangle$：多了一个**沿梯度方向的一阶项**。这不是修饰 —— 正因为多这一项，才能把 $\nabla F$ 从式子里消掉、使极值原理可用。**Bochner 论证的固定形状就是：二阶算子 + 一阶修正项 + 关于 $\lvert\nabla^2\rvert^2$ 的代数不等式**，缺一不可。这也是 [[04-调和函数·Harnack不等式与内部梯度估计]] §4.2 与 §24 篇末提到的"Bochner 技巧"的真身。

> **补充（无界域上的收尾要小心）**　第四步默认了「$F$ 的最大值在内部取到」。在**有界**区域上这没问题（配合边界处理）；但在 $\mathbb R^n$ 上，$\sup F$ 未必取到，严格做法是：
>
> 1. 先用**局部**估计（§25.3）在 $B_R\times(0,T]$ 上得到带 $C/R^2$ 的界；
> 2. 令 $R\to\infty$，$C/R^2\to0$，得 (25.7)。
>
> 讲稿直接从"同样有"跳到结论，跳过的正是这一步。若额外假设 $\sup F<+\infty$ 并配合截断函数，也可以直接走极值原理，但那时要处理截断带来的 $\nabla\eta$ 项，计算反而更长。

---

## 25.3 局部梯度估计

**定理 25.3.1**　设 $R>0$，$u$ 是 $B_R\times(0,T]$ 上热方程的**正**解，$\alpha>1$。则对任意 $(x,t)\in B_{R/2}\times(0,T]$，

$$
\frac{\left|\nabla u\right|^2}{u^2}-\alpha\frac{\partial_tu}{u}\ \le\ \frac{n\alpha^2}{2t}+\frac{C\alpha^2}{R^2}\left(1+\frac{\alpha^2}{\alpha-1}\right),\tag{25.11}
$$

其中 $C$ 只依赖于 $n$。

讲稿在此写「证明：略，用类似 Laplace 方程的 Bochner 技巧」。**这是全篇唯一整节跳过证明的地方**，而它偏偏是 §25.4 Harnack 的**唯一入口**。所以至少要说清它的三个部件，否则读者拿不到任何可复现的信息：

| 部件 | 作用 | 与全局情形的差别 |
|---|---|---|
| Bochner 计算 | 得到 (25.10) 那样的演化不等式 | 完全相同 |
| **截断函数** $\eta\in C_c^\infty(B_R)$、$\eta\equiv1$ 于 $B_{R/2}$ | 把 $\eta^2F$ 当待估量，让极值原理在**有界**区域上可用 | 全局情形不需要 |
| $\Delta\eta$、$\lvert\nabla\eta\rvert$ 的估计 | 截断产生的项被 $\frac{C}{R^2}$ 吸收，同时吃掉 $\lvert\nabla^2f\rvert^2$ 的一部分 | 产生 $C/R^2$ 这一项 |

三个参数各司其职：$R$ 决定"局部化到什么尺度"（$C/R^2$ 就是这一层的代价）；$\alpha$ 是可调权重，$\alpha>1$ 是计算中出现的门槛（$\frac{\alpha^2}{\alpha-1}$ 在 $\alpha\downarrow1$ 时发散）；$\frac{n\alpha^2}{2t}$ 在 $\alpha\to0$ 时趋于 $0$，与 $\alpha=0$ 的"纯梯度估计"吻合。

> **要点（$\alpha=2$ 不是随便取的）**　§25.4 的 Harnack 证明要对 $\log u$ 用 $\lvert\nabla f\rvert^2-2f_t$，即 $\alpha=2$。为什么必须是 $2$？因为在 Harnack 的积分里要把 $\lvert\nabla f\rvert\lvert x-y\rvert-\frac{1}{\alpha}\lvert\nabla f\rvert^2(t_2-t_1)$ 用 AM–GM 压成 $\frac{\alpha\lvert x-y\rvert^2}{4(t_2-t_1)}$，而**只有 $\alpha=2$ 时这个系数变成干净的 $\frac12$**，与后面 $\frac{1}{2(t_2-t_1)}$ 正好配成 $\frac{\lvert x-y\rvert^2}{2(t_2-t_1)}$。换成别的 $\alpha$ 只是常数变差，取 $2$ 是为了式子好看 —— 但这一步"为了好看"的选择在后面每一步都要用到，不能中途改。

---

## 25.4 Harnack 不等式

**定理 25.4.1（Harnack）**　设 $u$ 是 $B_R\times(0,T]$ 上热方程的**正**解。则对任意 $0<t_1<t_2\le T$ 与 $x,y\in B_{R/2}$，

$$
u(x,t_1)\ \le\ u(y,t_2)\left(\frac{t_2}{t_1}\right)^{n}\exp\left(\frac{R^2}{2(t_2-t_1)}+\frac{C\left(t_2-t_1\right)}{R^2}\right).\tag{25.12}
$$

> **理解（这个不等式为什么长这样）**　它说：**同一时刻、同一小球内**的两点值可比（[[04-调和函数·Harnack不等式与内部梯度估计]] 的椭圆 Harnack）；而抛物情形可以比较**不同时刻**的两点，代价是两个指数因子 ——
>
> - $\left(\frac{t_2}{t_1}\right)^n$：热核在 $t$ 时刻的峰值是 $(4\pi t)^{-n/2}$，随时间**衰减**，故"晚"的那一点天然更小，要补上 $(t_2/t_1)^{n/2}$ 的平方；
> - $e^{\frac{R^2}{2(t_2-t_1)}}$：相隔时间越短，两点之间能维持的梯度越陡，故只有 $\frac{1}{t_2-t_1}$ 型的控制；
> - $e^{\frac{C(t_2-t_1)}{R^2}}$：长时间累积的误差项。
>
> **方向不能反**：$t_1<t_2$ 时"早"的被"晚"的控制。反过来问"晚的能否被早的控制"，答案是**不能**（可以差任意倍，因为热核峰值随 $t$ 下降）。这是抛物 Harnack 与椭圆 Harnack 最本质的区别。

**证明**　记 $f:=\log u$。由定理 25.3.1 取 $\alpha=2$：

$$
\left|\nabla f\right|^2-2f_t\ \le\ \frac{2n}{t}+\frac{C}{R^2}\qquad\text{于 }B_{R/2}\times(0,T].\tag{25.13}
$$

（把 (25.11) 里的 $\frac{C\alpha^2}{R^2}\left(1+\frac{\alpha^2}{\alpha-1}\right)=\frac{20C}{R^2}$ 吸收进 $\frac{C}{R^2}$，$C$ 重新命名 —— 见订正第 4 条。）

**沿路径积分。** 取从 $(y,t_2)$ 到 $(x,t_1)$ 的直线路径

$$
x_s=y+s(x-y),\qquad t_s=t_2+s(t_1-t_2),\qquad 0\le s\le1
$$

（$t_s$ 从 $t_2$ 递减到 $t_1$；全程 $x_s\in B_{R/2}$，因为 $B_{R/2}$ 凸而 $x,y$ 都在其中）。由微积分基本定理，

$$
\log\frac{u(x,t_1)}{u(y,t_2)}=f(x,t_1)-f(y,t_2)=\int_0^1\frac{d}{ds}f(x_s,t_s)\,ds=\int_0^1\left[\left\langle\nabla f(x_s,t_s),\,x-y\right\rangle+\left(t_1-t_2\right)\partial_tf(x_s,t_s)\right]ds .
$$

记 $\tau:=t_2-t_1>0$，则 $t_1-t_2=-\tau$，逐点估计被积函数：

$$
\left\langle\nabla f,x-y\right\rangle-\tau f_t\ \le\ \left|\nabla f\right|\left|x-y\right|-\frac{\tau}{2}\left|\nabla f\right|^2+\tau\left(\frac{n}{t_s}+\frac{C}{2R^2}\right),
$$

其中用了 (25.13) 的等价形式 $-f_t\le\frac{n}{t_s}+\frac{C}{2R^2}-\frac12\lvert\nabla f\rvert^2$。

**AM–GM 压掉 $\lvert\nabla f\rvert$。** 对 $a=\lvert\nabla f\rvert$ 配方：

$$
\left|\nabla f\right|\left|x-y\right|-\frac{\tau}{2}\left|\nabla f\right|^2=\frac{\left|x-y\right|^2}{2\tau}-\frac{\tau}{2}\left(\left|\nabla f\right|-\frac{\left|x-y\right|}{\tau}\right)^2\ \le\ \frac{\left|x-y\right|^2}{2\tau}.
$$

**对 $s$ 积分。** 剩下两项：

$$
\int_0^1\frac{ds}{t_2+s(t_1-t_2)}=\frac{1}{t_1-t_2}\ln\frac{t_1}{t_2}=\frac{1}{\tau}\ln\frac{t_2}{t_1}\quad\Longrightarrow\quad \tau n\int_0^1\frac{ds}{t_s}=n\ln\frac{t_2}{t_1} .
$$

合并：

$$
\log\frac{u(x,t_1)}{u(y,t_2)}\le\frac{\left|x-y\right|^2}{2\tau}+n\ln\frac{t_2}{t_1}+\frac{C\tau}{2R^2}\le\frac{R^2}{2\tau}+n\ln\frac{t_2}{t_1}+\frac{C\tau}{2R^2},
$$

末步用 $\lvert x-y\rvert\le R$（$x,y\in B_{R/2}$）。指数化即 (25.12)。$\square$

---

## 25.5 强极值原理

**定理 25.5.1（强极值原理）**　设 $u$ 是 $\Omega\times(0,T]$ 上热方程的解。若存在 $(x_0,t_0)\in\Omega\times(0,T]$ 使

$$
u(x_0,t_0)=\max_{\overline{\Omega\times(0,T]}}u ,
$$

则 $u$ 在 $\Omega\times(0,t_0]$ 上恒为常数。

讲稿到此为止，**没有证明**。这一步其实就是 Harnack 不等式的直接应用 —— 而这也解释了它为什么被排在 §25.4 之后（尽管它属于「极值原理」这一族结论，归 [[24-弱极值原理与热方程初边值问题的唯一性]] §24.6 那张表）。

**证明（补充）**　记 $M:=u(x_0,t_0)$、$w:=M-u\ge0$。$w$ 同样解热方程，$w(x_0,t_0)=0$。要证 $w\equiv0$ 于 $\Omega\times(0,t_0]$。

**第一步：在一个小球柱上全零。** 取 $R>0$ 使 $B_R(x_0)\subseteq\Omega$（$x_0\in\Omega$ 开）。

对任意 $\varepsilon>0$，$w+\varepsilon\ge\varepsilon>0$ 于 $B_R\times(0,T]$ 且解热方程，故可用定理 25.4.1。任取 $t\in(0,t_0)$ 与 $x\in B_{R/2}(x_0)$，以 $(x,t)$ 为"早"、$(x_0,t_0)$ 为"晚"：

$$
w(x,t)+\varepsilon\ \le\ \left(w(x_0,t_0)+\varepsilon\right)\left(\frac{t_0}{t}\right)^{n}e^{\frac{R^2}{2(t_0-t)}+\frac{C(t_0-t)}{R^2}}=\varepsilon\cdot C_{t,t_0,R} .
$$

（第一个等号用到 $w(x_0,t_0)=0$。常数 $C_{t,t_0,R}$ 对**固定的** $t<t_0$ 有限 —— 这正是必须取 $t<t_0$ 而不能取 $t=t_0$ 的原因。）令 $\varepsilon\to0^+$：

$$
w(x,t)=0\quad\Longrightarrow\quad u(x,t)=M,\qquad \forall x\in B_{R/2}(x_0),\ \forall t\in(0,t_0).\tag{25.14}
$$

**第二步：沿空间传播。** 令

$$
A:=\left\{x\in\Omega:\ u(x,t)=M\ \text{对一切 }t\in(0,t_0]\right\}.
$$

- $A\ne\varnothing$：由 (25.14) 与 $u$ 的连续性，$x_0\in A$；
- $A$ 在 $\Omega$ 中**开**：对 $x_1\in A$，取 $R_1>0$ 使 $B_{R_1}(x_1)\subseteq\Omega$，把第一步的论证（以 $x_1$ 为中心、$(x_1,t_0)$ 为"晚"点）重做一遍 —— 注意 $u(x_1,t_0)=M$ 是最大值，论证原封不动 —— 得 $B_{R_1/2}(x_1)\subseteq A$；
- $A$ 在 $\Omega$ 中**闭**：$u$ 连续，$A$ 是闭集的原像之交。

$\Omega$ 连通且 $A$ 非空、既开又闭，故 $A=\Omega$。$\square$

> **要点（为什么只能推到 $t_0$）**　Harnack 只会**往前**传（$t_1<t_2$），所以从 $(x_0,t_0)$ 出发只能回溯到 $t<t_0$，永远够不到 $t>t_0$。这不是证明的缺陷，而是事实：$t>t_0$ 的值还没被方程决定，完全可以掉下来。**对照椭圆情形**：Laplace 的强极值原理说"内部取到极值 ⟹ 全 $\Omega$ 常数"，那里没有时间，所以没有这个限制。这是 [[24-弱极值原理与热方程初边值问题的唯一性]] §24.6 那张表里唯一"抛物更弱"的一格。

---

## 25.6 Harnack 的推论：区域内部的整体比较

**定理 25.6.1**　设 $\Omega\subseteq\mathbb R^n$ 为有界光滑区域，$u>0$ 于 $\Omega_T$ 且解热方程。则对任意 $\Omega'\Subset\Omega$ 与 $0<t_1<t_2\le T$，

$$
\sup_{x\in\Omega'}u(x,t_1)\ \le\ C\inf_{x\in\Omega'}u(x,t_2),\qquad C=C\left(n,\Omega',\Omega,t_1,t_2\right).\tag{25.15}
$$

定理 25.4.1 是**局部**的：一个球、两个点、常数依赖球半径 $R$。推论把它升级成**整体**的：整个内层区域、两个时刻、常数只依赖区域的几何与这两个时刻。这一步是 §25.7 强极值原理的唯一入口 —— 强极值原理要在"任取的 $\Omega'\Subset\Omega$"上做，而定理 25.4.1 给不出 $\Omega'$ 上的 $\sup\le C\inf$，只能给出球上的。

**证明**　四步。

**第一步：退到内层。** Harnack 要求整球 $B_R$ 落在 $\Omega$ 里，而 $\Omega'$ 的边界点没有这个余量。故先退一层：记

$$
d(y,\partial\Omega):=\inf_{z\in\partial\Omega}\lvert y-z\rvert,\qquad \Omega_\delta:=\left\{y\in\Omega:\ d(y,\partial\Omega)>\delta\right\}.
$$

$\Omega$ 光滑保证 $\delta>0$ 充分小时 $\Omega_\delta$ 仍连通，且 $\Omega'\subset\Omega_\delta$。取 $x_1,x_2\in\overline{\Omega_\delta}$ 使

$$
u(x_1,t_1)=\sup_{x\in\overline{\Omega_\delta}}u(x,t_1),\qquad u(x_2,t_2)=\inf_{x\in\overline{\Omega_\delta}}u(x,t_2).
$$

只需证 $u(x_1,t_1)\le Cu(x_2,t_2)$：因为 $\sup_{\Omega'}\le\sup_{\Omega_\delta}$、$\inf_{\Omega_\delta}\le\inf_{\Omega'}$，中间这一步成立则 (25.15) 自动成立。而且证明实际给出的比 (25.15) 更强 —— $\inf$ 是取在**更大的** $\overline{\Omega_\delta}$ 上的。

**第二步：空间上走球链。** 取 $0<R<\delta$。$\overline{\Omega_\delta}$ 紧且连通，故存在有限个点 $y_0,\dots,y_N\in\overline{\Omega_\delta}$，使 $y_0=x_1$、$y_N=x_2$，且

$$
\overline{\Omega_\delta}\subset\bigcup_{i=0}^{N}B_{R/2}(y_i),\qquad \lvert y_i-y_{i-1}\rvert<\frac R2,\qquad i=1,\dots,N .
$$

这条链要同时满足两件事：$B_R(y_i)\subset\Omega$（因 $d(y_i,\partial\Omega)\ge\delta>R$，这样才能用定理 25.4.1），以及 $y_{i-1},y_i\in B_{R/2}(y_i)$（因 $\lvert y_i-y_{i-1}\rvert<R/2$，这样两端都落在 Harnack 允许的"小半圈"里）。

> **要点（$N$ 必须与 $x_1,x_2$ 无关）**　做法是先固定一个半径 $R/4$ 的有限球覆盖（$\overline{\Omega_\delta}$ 紧），再沿相交的球挑出链上的点 —— 于是 $N$ 只依赖 $n,\Omega_\delta,\delta,R$。这不是技术性附注：下面第四步算出的因子是 $\exp\left(\frac{N^2R^2}{2(t_2-t_1)}\right)$，若 $N$ 随别的东西（比如时间的细分）无限增大，这个因子会炸掉，整个推论作废。

**第三步：时间上等分。** 令

$$
s_i=t_1+\frac iN\left(t_2-t_1\right),\qquad i=0,\dots,N .
$$

以 $B_R(y_i)$ 为基准球、$(y_{i-1},s_{i-1})$ 为"早"、$(y_i,s_i)$ 为"晚"，定理 25.4.1 给出

$$
u(y_{i-1},s_{i-1})\le C_i\,u(y_i,s_i),\qquad C_i=\left(\frac{s_i}{s_{i-1}}\right)^{n}\exp\left(\frac{R^2}{2(s_i-s_{i-1})}+C_n\frac{s_i-s_{i-1}}{R^2}\right).
$$

代入 $s_i-s_{i-1}=(t_2-t_1)/N$：

$$
C_i=\left(\frac{s_i}{s_{i-1}}\right)^{n}\exp\left(\frac{NR^2}{2(t_2-t_1)}+C_n\frac{t_2-t_1}{NR^2}\right).\tag{25.16}
$$

**第四步：常数核算。** 讲稿写到 $\prod_iC_i$ 就收了，这里算完。三段因子分开处理：

(a) **时间比。** 由 $s_{i-1}\ge t_1$ 得 $\dfrac{s_i}{s_{i-1}}=1+\dfrac{(t_2-t_1)/N}{s_{i-1}}\le1+\dfrac{t_2-t_1}{Nt_1}$，故

$$
\prod_{i=1}^{N}\left(\frac{s_i}{s_{i-1}}\right)^{n}\le\left(1+\frac{t_2-t_1}{Nt_1}\right)^{nN}\le\exp\left(\frac{n(t_2-t_1)}{t_1}\right),
$$

末步用 $(1+a/N)^N\le e^a$。

(b) **主指数。** $N$ 步完全相同，乘起来

$$
\prod_{i=1}^{N}\exp\left(\frac{NR^2}{2(t_2-t_1)}\right)=\exp\left(\frac{N^2R^2}{2(t_2-t_1)}\right).
$$

注意这里是 $N^2$ 不是 $N$：单步的指数是 $\frac{R^2}{2\Delta s}$，而 $\Delta s=\frac{t_2-t_1}{N}$ 已经把 $N$ 送进了分母，再累乘 $N$ 步就平方了。

(c) **误差项。** $\displaystyle\prod_{i=1}^{N}\exp\left(C_n\frac{t_2-t_1}{NR^2}\right)=\exp\left(C_n\frac{t_2-t_1}{R^2}\right)$，$N$ 恰好抵消。

合并三段：

$$
u(x_1,t_1)\le\exp\left(\frac{n(t_2-t_1)}{t_1}+\frac{N^2R^2}{2(t_2-t_1)}+\frac{C_n(t_2-t_1)}{R^2}\right)u(x_2,t_2).
$$

右端唯一含 $R$ 的地方就是指数里的两项，而 $R$ 在第二步已固定、$N=N(R,\Omega_\delta)$ 随之固定，于是整个指数只依赖 $n,\Omega',\Omega,t_1,t_2$。取它作 $C$ 即得 (25.15)。$\square$

> **要点（为什么"空间退一层 + 时间分 $N$ 份"必须同时做）**　Harnack 一次只能比"同一球内、两个**严格不同**时刻"的两个点。要从 $(x_1,t_1)$ 走到 $(x_2,t_2)$，空间要走 $N$ 段、时间也要走 $N$ 段，两者必须**同步**：第 $i$ 步同时把位置从 $y_{i-1}$ 推到 $y_i$、把时间从 $s_{i-1}$ 推到 $s_i$。若只分空间不分时间（$s_i\equiv t_2$），第一步的"晚"就是终点，后面几步无时间可走；若只分时间不分空间，又跨不过区域。**抛物 Harnack 的路径是"斜着走"的** —— 这与椭圆 Harnack 只需在空间里连一条折线，是最本质的技术差别。

---

## 25.7 强极值原理：讲稿的另一种证法

§25.5 走的是「$w=M-u$ + 小球柱 + 既开又闭」。讲稿走另一条路：**把 $u$ 直接补成正解，然后一次用 §25.6 的推论**。两条路的实质相同，但讲稿版更短，也更清楚地点出了「强极值原理 = Harnack + $\varepsilon$ 扰动」这句话。

**定理 25.7.1**（同定理 25.5.1）　设 $u$ 在 $\overline{\Omega_T}$ 上连续、在 $\Omega_T$ 内解热方程。若存在 $(x_0,t_0)\in\Omega_T$ 使 $u(x_0,t_0)=\sup_{\Omega_T}u$，则 $u$ 在 $\Omega_{t_0}=\Omega\times(0,t_0]$ 上恒为常数。

**证明（讲稿版）**　对任意 $\varepsilon>0$ 令

$$
v_\varepsilon(x,t):=\sup_{\Omega_T}u-u(x,t)+\varepsilon>0 .
$$

$v_\varepsilon$ 同样解热方程，且在最大值点处 $v_\varepsilon(x_0,t_0)=\varepsilon$。

任取含 $x_0$ 的子区域 $\Omega'\Subset\Omega$ 与 $0<t_1<t_0$。对 $v_\varepsilon$ 用定理 25.6.1（区域 $\Omega'$，时刻 $t_1<t_0$）：

$$
\sup_{x\in\Omega'}v_\varepsilon(x,t_1)\ \le\ C\inf_{x\in\Omega'}v_\varepsilon(x,t_0)\ \le\ C\,v_\varepsilon(x_0,t_0)=C\varepsilon ,
$$

中间一步用了 $x_0\in\Omega'$。把 $v_\varepsilon$ 的定义代回去，即对一切 $x\in\Omega'$，

$$
u(x,t_1)\ \ge\ \sup_{\Omega_T}u-(C-1)\varepsilon .
$$

令 $\varepsilon\to0^+$ 得 $u(x,t_1)\equiv\sup_{\Omega_T}u$ 于 $\Omega'$。由 $\Omega'\Subset\Omega$ 与 $t_1\in(0,t_0)$ 的任意性，再配合 $u$ 的连续性，得 $u\equiv\sup_{\Omega_T}u$ 于 $\Omega_{t_0}$。$\square$

> **要点（三处容易漏掉的细节）**
>
> 1. **必须 $t_1<t_0$ 严格。** 定理 25.6.1 的常数含 $\frac{1}{t_2-t_1}$，取 $t_2=t_0$ 而让 $t_1\to t_0$ 就发散。这正是"强极值原理只能推到 $t_0$、推不到 $t>t_0$"的技术根源 —— 不是证明偷懒，是 Harnack 本身有方向。
> 2. **$\varepsilon$ 不能省。** $\sup u-u$ 只是**非负**，在 $u$ 取到最大值的点为零，不能作除数、Harnack 的"正解"前提不成立。$+\varepsilon$ 是把非负抬成正解的全部代价，最后再让 $\varepsilon\to0$ 收回。
> 3. **这里用的是 §25.6 的推论，不是 §25.4 的定理。** 定理 25.4.1 只在**一个球**内比较两点，而 $\Omega'$ 一般不是球。讲稿直接从"由 Harnack 不等式"跳到 $\Omega'$ 上的 $\sup\le C\inf$，中间那一步（球链）正是 §25.6 补的。

> **对照（两条路的分工）**　§25.5 的"既开又闭"是在**空间**上做延拓：先把一片区域证成常数，再靠连通性铺满整个 $\Omega$，时间方向上是逐点固定的。§25.7 则是直接在**整个 $\Omega'$ 上**用一次整体比较，把所有空间点一次性拿下。前者只需要球版 Harnack，后者需要球链版 —— **省下来的篇幅，代价是背后那条球链**。两条路都证明了同一件事，选哪条只看手头有哪个版本的 Harnack。

---

## 25.8 与调和函数的对照

| | 调和函数（[[04-调和函数·Harnack不等式与内部梯度估计]]） | 热方程（本篇） |
|---|---|---|
| 梯度估计 | $\lvert\nabla\log u\rvert\le C/R$ | $\lvert\nabla\log u\rvert^2-\partial_t\log u\le n/2t$ |
| 估计的"尺度参数" | $R$（空间） | $t$（时间） |
| 证明路线 A | Poisson 核 + Cauchy–Schwarz | 热核 + 加权 Cauchy–Schwarz（§25.1） |
| 证明路线 B | Bochner + 极值原理 | Bochner + 极值原理（§25.2） |
| 常数的来源 | 核的 $\lvert\nabla H\rvert/H\sim1/R$ | 核自身取等号（引理 25.1.2） |
| Harnack | $u(x)\le Cu(y)$，**同层** | $u(x,t_1)\le Cu(y,t_2)$，**只能早 $\le$ 晚**（§25.4） |
| Harnack 的常数 | 只依赖 $n$ 与 $R$ | 多出 $\left(\frac{t_2}{t_1}\right)^n$ 与 $\frac{R^2}{t_2-t_1}$ |
| 强极值原理 | 内部极值 ⟹ 全区域常数 | 内部极值 ⟹ $\Omega\times(0,t_0]$ 上常数（§25.5） |

**一句话总结这张表**：抛物情形把椭圆的"空间尺度 $R$"换成"时间尺度 $t$"，并在每个结论里多一个 $\partial_t$ 或一个时间因子；方向性（只能往前）是唯一的定性差别。

---

## 订正

1. **局部梯度估计的区域写成 $B_2$，应为 $B_{R/2}$。** (25.11) 右端含 $\frac{C}{R^2}$，说明结论是在"比定义域小一圈"的球上成立的；而 §25.4 的 Harnack 证明里路径 $x_s$ 必须全程落在这个球里（因为要用 (25.13)），所以它必须是 $B_{R/2}$。$B_2$ 这个写法与满篇的 $R$ 无法对应 —— 若照字面读成"半径 $2$ 的球"，则 $R$ 与区域无关，$\frac{C}{R^2}$ 就失去了"局部化代价"的含义。判据：**$R$ 必须同时出现在区域名与界里**，且区域尺度是 $R$ 的常数倍。

2. **「若 $u$ 为一般的正解，同样有梯度估计」这句没证，且直接对 $F$ 用极值原理在 $\mathbb R^n$ 上不严谨。** 见 §25.2 的补充：$\sup F$ 未必取到，严格做法是由局部估计令 $R\to\infty$。讲稿把一个需要"局部化 + 取极限"的结论压成一句话，读者若照着 §25.1 的表示法去证，会发现没有初值就没有那个概率测度 $\nu$，路走不通。本篇 §25.2 给了 Bochner 论证（可证），并注明两种收尾方式。

3. **全局梯度估计的定理条件里"连续"是多余的。** 定理 25.1.1 的证明只用到：$\varphi\ge0$、$\varphi$ 有界、$\varphi>0$ 于某个球（从而 $u>0$）。所以 $\varphi\in L^\infty$ 且 $\varphi>0$ 于某正测集就够了。另外"取 $y_0$ 使 $\varphi(y_0)>0$"这一句对可测函数不成立（点值无意义），应改为"存在 $r>0$ 与 $c>0$ 使 $\varphi\ge c$ 于某球"。$H_t*\varphi$ 的光滑性由 $H_t\in\mathcal S$ 保证，与 $\varphi$ 的连续性无关（[[22-Cauchy问题与热核]] 命题 22.6.1）。

4. **Harnack 证明里的 $C$ 与局部估计里的 $C$ 不是同一个数。** 由 (25.11) 取 $\alpha=2$ 得 $\frac{C\alpha^2}{R^2}\left(1+\frac{\alpha^2}{\alpha-1}\right)=\frac{20C}{R^2}$，讲稿直接写成 $\frac{C}{R^2}$，即 $C$ 被悄悄乘了 $20$。结论不受影响（常数无关紧要），但**若不说明就直接对照两处 $C$，会误以为可以逐项相等**。本篇 (25.13) 后已注明。

5. **强极值原理只给陈述、没有证明。** 它是 Harnack 的直接推论（§25.5），证明不到十行，省掉之后「为什么强极值原理排在 Harnack 之后」这个编排上的用意也一并丢失了。本篇 §25.5 已补，并点明「只能往前推到 $t_0$」的原因。

6. **记号不一致：定理 25.1.1 的陈述写 $u(t,x)=H(\cdot,t)*\varphi$，其后所有推导用 $u(x,t)$。** 热方程研究里 $(x,t)$ 是标准次序（空间在前），$u(t,x)$ 是概率论的写法（时间在前）。全篇应统一为 $u(x,t)$ —— 尤其因为 Harnack 里 $x$ 与 $t$ 的角色完全不同（$t$ 有方向性、$x$ 没有），次序写混会导致把"只能往前"记成"只能往外"。

7. **Harnack 推论的证明收在 $\prod_iC_i$ 上，常数没算完，也就没交代 $N$ 与 $R$ 的关系。** 由 (25.16)，单步因子含 $\exp\left(\frac{NR^2}{2(t_2-t_1)}\right)$，累乘 $N$ 步得 $\exp\left(\frac{N^2R^2}{2(t_2-t_1)}\right)$ —— 只要 $N$ 有一点点依赖 $R$，这个量就可能随参数发散。结论能成立，全靠第二步里"$N$ 由 $R/4$ 球覆盖确定、与 $x_1,x_2$ 无关"这句话。讲稿把它写成一句夹注（"故 $N$ 的上界可与 $x_1,x_2$ 无关"），但它其实是整条论证的命门。本篇 §25.6 第四步已补完，并把这个因子单独点出。

8. **推论的陈述把 $\sup/\inf$ 取在 $\Omega'$ 上，证明里却取在 $\overline{\Omega_\delta}$ 上，两处不一致。** 陈述用 $\Omega'$ 是对的（结论只该依赖 $\Omega'$），但证明真正得到的是 $\sup_{\Omega'}u(\cdot,t_1)\le C\inf_{\overline{\Omega_\delta}}u(\cdot,t_2)$ —— $\inf$ 的自然位置在更大的 $\overline{\Omega_\delta}$ 上，然后靠 $\inf_{\Omega_\delta}\le\inf_{\Omega'}$ 收到 $\Omega'$。讲稿没有写这一步"收"，读者会以为证明证的正是陈述本身。另外"$\Omega_\delta$ 连通"是借"$\Omega$ 光滑、$\delta$ 充分小"得到的，讲稿写对了，但没说明为什么需要连通 —— 因为球链必须能在 $\Omega_\delta$ 里连起来，不连通时 $x_1,x_2$ 可能落在不同分支上。

9. **§3 强极值原理的证明写"由 Harnack 不等式"，但球版 Harnack 直接给不出 $\Omega'$ 上的 $\sup\le C\inf$。** 实际用的是 §2 的推论（球链版）。严格说 §3 依赖 §2，而 §2 的推论依赖 §2 的定理 —— 这个层层依赖讲稿没标出来，照字面读会以为一步就到位。本篇 §25.7 已把依赖关系写明，并给出与 §25.5 的"$w=M-u$ + 既开又闭"证法的对照。

---

系列：[[24-弱极值原理与热方程初边值问题的唯一性]] · [[26-Cauchy问题的唯一性与Tychonoff反例]] · [[04-调和函数·Harnack不等式与内部梯度估计]] · [[22-Cauchy问题与热核]]
