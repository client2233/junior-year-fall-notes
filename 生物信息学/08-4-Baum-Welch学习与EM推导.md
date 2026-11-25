# 08-4 Baum-Welch 学习与 EM 推导

> 系列：[[01-生物学基础]] · [[02-生物序列数据库与检索]] · [[03-1-全局比对·Needleman-Wunsch]] · [[03-2-局部比对与仿射罚分]] · [[03-3-启发式算法·FASTA与BLAST]] · [[04-速查·矩阵与公式]] · [[05-1-Linux命令行与测序数据]] · [[05-2-Python序列处理·从DNA到蛋白质]] · [[05-3-结构与组学数据·PDB与单细胞]] · [[06-1-机器学习概念·PCA·聚类与降维]] · [[06-2-监督学习·损失函数·梯度下降与数据集划分]] · [[06-3-分类评价指标·混淆矩阵·ROC与PR]] · [[06-4-树模型与ALL-AML实战·数据泄漏]] · [[06-5-感知机·原始形式与对偶形式]] · [[07-1-分子进化·系统树与距离矩阵]] · [[07-2-UPGMA·超度量与分子钟]] · [[07-3-Neighbor-Joining·最小演化推导与算例]] · [[07-4-多序列比对与ClustalW]] · [[07-5-PSI-BLAST·PSSM与位置特异性打分]] · [[07-6-sequence-profile与profile比对]] · [[08-1-隐马尔可夫模型·定义与三个基本问题]] · [[08-2-概率计算·前向与后向算法]] · [[08-3-Viterbi解码·动态规划求最优路径]] · [[08-5-profile-HMM与HMMER]]

问题②：给定观测序列 $O$（可能还有状态序列 $I$），估计 $\lambda=(A,B,\pi)$ 使 $P(O\mid\lambda)$ 最大。

两种情形，难度差一个量级：

| 情形 | 已知 | 方法 | 有没有闭式解 |
|---|---|---|---|
| 一（监督） | $O$ 和 $I$ | 直接数频数 | 有，一步到位 |
| 二（无监督） | 只有 $O$ | Baum–Welch（EM） | 没有，迭代 |

## 24.1 两个后验量：$\gamma$ 与 $\xi$

Baum–Welch 的全部内容就是「用当前的 $\lambda'$ 算出后验期望计数，再拿计数当新参数」。两个量：

**$\gamma_t(i)$**：给定观测，时刻 $t$ 处于 $S_i$ 的概率

$$
\gamma_t(i)=P\bigl(q_t=S_i\mid O,\lambda\bigr)=\frac{P(O,q_t=S_i\mid\lambda)}{P(O\mid\lambda)}
$$

用 $\alpha,\beta$ 表示（这是唯一实用的算法）：

$$
\gamma_t(i)=\frac{\alpha_t(i)\,\beta_t(i)}{P(O\mid\lambda)}
=\frac{\alpha_t(i)\,\beta_t(i)}{\sum_{i=1}^{N}\alpha_t(i)\,\beta_t(i)}
$$

分母那样写是为了**免去单独算 $P(O\mid\lambda)$**，同时天然归一（$\sum_i\gamma_t(i)=1$）。

**$\xi_t(i,j)$**：时刻 $t$ 在 $S_i$、$t+1$ 在 $S_j$ 的概率

$$
\xi_t(i,j)=P\bigl(q_t=S_i,\,q_{t+1}=S_j\mid O,\lambda\bigr)
=\frac{\alpha_t(i)\,a_{ij}\,b_j(o_{t+1})\,\beta_{t+1}(j)}{P(O\mid\lambda)}
$$

分子四段的含义：**走到 $i$（$\alpha$）→ 跳到 $j$（$a_{ij}$）→ 在 $j$ 发射 $o_{t+1}$（$b$）→ 从 $j$ 看完剩下（$\beta$）**。这四段正好覆盖整条序列，不重不漏 —— 看图记比背公式可靠。

归一：$\sum_i\sum_j\xi_t(i,j)=1$（对每个 $t$ 都是概率分布）。

**两者的关系**

$$
\sum_{j=1}^{N}\xi_t(i,j)=\gamma_t(i),\qquad t=1,\dots,T-1
$$

> **订正** 课件第 39 页把这条写成 $\gamma_t(i)=\sum_{j=1}^{N}\xi_t(i,j)$ 而**没有限定 $t\le T-1$**。$t=T$ 时 $\xi_T$ 根本没有定义（没有 $t+1$），此时只能 $\gamma_T(i)=\alpha_T(i)/P(O\mid\lambda)$。写代码时若无条件套这个式子，最后一列 $\gamma$ 会全变成 0。

## 24.2 期望计数

| 量 | 表达式 | 读法 |
|---|---|---|
| 状态 $S_i$ 出现的期望次数 | $\sum_{t=1}^{T}\gamma_t(i)$ | 整条序列里在 $i$ 待了几步 |
| 从 $S_i$ **出发**的转移期望次数 | $\sum_{t=1}^{T-1}\gamma_t(i)$ | 少了最后一步（$T$ 时无出发） |
| $S_i\to S_j$ 的转移期望次数 | $\sum_{t=1}^{T-1}\xi_t(i,j)$ | 逐边计数 |

第二行与第三行必须一致：$\sum_j\sum_t\xi_t(i,j)=\sum_t\gamma_t(i)$（$t\le T-1$）。这个恒等式是 Baum–Welch 分母的合法性来源，值得记。

> **补充** 课件第 40 页把第二行叫「状态 $S_i$ 转移的期望数」，措辞含糊（听起来像「转移到 $S_i$ 的次数」）。准确说法是「**从** $S_i$ 出发的转移次数」。中文里「A 转移」有歧义，英文 departures from $S_i$ 没有。

## 24.3 情形一：已知 $O$ 和 $I$（监督学习）

这是最简单的情形，纯数频数：

$$
\hat\pi_i=\mathbb 1\{i_1=S_i\}
$$

$$
\hat a_{ij}=\frac{A_{ij}}{\sum_{j=1}^{N}A_{ij}},\qquad A_{ij}=\#\{t: i_t=S_i,\ i_{t+1}=S_j\}
$$

$$
\hat b_j(k)=\frac{B_{jk}}{\sum_{k=1}^{M}B_{jk}},\qquad B_{jk}=\#\{t: i_t=S_j,\ o_t=v_k\}
$$

即「$i\to j$ 的次数 / 从 $i$ 出发的总次数」与「$j$ 发射 $k$ 的次数 / 在 $j$ 的总次数」。

**坑**：分子为 0 时（某个转移从未出现），$\hat a_{ij}=0$。这在监督学习里是合法的（数据说它不可能），但会让后续前向概率整块为 0；若 $N$ 大、$T$ 小，大量 0 会把模型锁死。实践上加**平滑**：

$$
\hat a_{ij}=\frac{A_{ij}+\alpha}{\sum_j A_{ij}+N\alpha}
$$

$\alpha=1$ 就是拉普拉斯平滑。课件没提，但这是文本/序列标注的标准做法。

## 24.4 情形二：只有 $O$ —— 重估公式

把 24.3 里的「次数」换成「期望次数」：

$$
\boxed{\ \hat\pi_i=\gamma_1(i)\ }
$$

$$
\boxed{\ \hat a_{ij}=\frac{\sum_{t=1}^{T-1}\xi_t(i,j)}{\sum_{t=1}^{T-1}\gamma_t(i)}\ }
$$

$$
\boxed{\ \hat b_j(k)=\frac{\sum_{t=1,\,o_t=v_k}^{T}\gamma_t(j)}{\sum_{t=1}^{T}\gamma_t(j)}\ }
$$

三个式子的读法完全一致：**分子 = 「这件事发生了多少次」的期望，分母 = 「这个状态出现了多少次」的期望**。

- $\hat a_{ij}$ 的分母是 $\sum_{t=1}^{T-1}\gamma_t(i)$，**不含 $t=T$** —— 因为 $T$ 时刻不出发；
- $\hat b_j(k)$ 的分子只累加「观测恰为 $v_k$」的那些 $t$，分母却是全部 $T$；
- 三者都自动满足归一（分子之和 = 分母），不需要再归一化。

## 24.5 为什么是这三条：EM 推导

上面的公式是「猜」出来的（用期望代替计数）。严格证明要走一遍 EM，这条推导也是全讲唯一需要动笔的部分。

### 24.5.1 把 HMM 看成含隐变量的模型

$O$ 是观测，$I$（状态序列）是隐变量。完全数据的对数似然

$$
\log P(O,I\mid\lambda)=\log\pi_{i_1}+\sum_{t=2}^{T}\log a_{i_{t-1}i_t}+\sum_{t=1}^{T}\log b_{i_t}(o_t)
$$

**这个分解是全部推导的地基**：三个参数 $\pi,A,B$ **分别**只出现在一项里，所以后面能拆成三个独立的优化问题。

### 24.5.2 E 步：构造 $Q$ 函数

$$
Q(\lambda,\lambda')=\sum_{I}P(I\mid O,\lambda')\log P(O,I\mid\lambda)
$$

$\lambda'$ 是当前参数（固定），$\lambda$ 是待优化的。因为 $P(I\mid O,\lambda')=P(O,I\mid\lambda')/P(O\mid\lambda')$ 且分母与 $\lambda$ 无关，等价于优化

$$
Q(\lambda,\lambda')\ \sim\ \sum_{I}P(O,I\mid\lambda')\log P(O,I\mid\lambda)
$$

把 24.5.1 的分解代入，$Q$ 裂成三项：

$$
Q=\underbrace{\sum_{I}\log\pi_{i_1}P(O,I\mid\lambda')}_{\text{只含 }\pi}
+\underbrace{\sum_{I}\Bigl(\sum_{t=1}^{T-1}\log a_{i_ti_{t+1}}\Bigr)P(O,I\mid\lambda')}_{\text{只含 }A}
+\underbrace{\sum_{I}\Bigl(\sum_{t=1}^{T}\log b_{i_t}(o_t)\Bigr)P(O,I\mid\lambda')}_{\text{只含 }B}
$$

三块可以**分别最大化** —— 这就是 24.5.1 那句「分别出现在三项里」的全部价值。

### 24.5.3 M 步（一）：$\pi$

第一项里只有 $i_1$ 有关，把对 $I$ 的求和压到 $i_1$ 上：

$$
\sum_{I}\log\pi_{i_1}P(O,I\mid\lambda')=\sum_{i=1}^{N}\log\pi_i\,P(O,q_1=S_i\mid\lambda')
$$

约束 $\sum_i\pi_i=1$，拉格朗日函数

$$
L(\pi,\eta)=\sum_{i=1}^{N}\log\pi_i\,P(O,q_1=S_i\mid\lambda')+\eta\Bigl(\sum_{i=1}^{N}\pi_i-1\Bigr)
$$

求导置零：

$$
\frac{\partial L}{\partial\pi_i}=\frac{1}{\pi_i}P(O,q_1=S_i\mid\lambda')+\eta=0
\ \Longrightarrow\
P(O,q_1=S_i\mid\lambda')+\pi_i\eta=0
$$

两边对 $i$ 求和，注意 $\sum_i\pi_i=1$：

$$
P(O\mid\lambda')+\eta=0\ \Longrightarrow\ \eta=-P(O\mid\lambda')
$$

代回：

$$
\hat\pi_i=\frac{P(O,q_1=S_i\mid\lambda')}{P(O\mid\lambda')}=\gamma_1(i)
$$

$\square$

### 24.5.4 M 步（二）：$A$

第二项把求和压到 $(i_t,i_{t+1})$ 上：

$$
\sum_{i=1}^{N}\sum_{j=1}^{N}\sum_{t=1}^{T-1}\log a_{ij}\,P(O,q_t=S_i,q_{t+1}=S_j\mid\lambda')
$$

约束是**逐行**的（$\sum_j a_{ij}=1$，$N$ 个约束，$N$ 个乘子 $\eta_i$）。对 $a_{ij}$ 求导：

$$
\frac{1}{a_{ij}}\sum_{t=1}^{T-1}P(O,q_t=S_i,q_{t+1}=S_j\mid\lambda')+\eta_i=0
$$

两边对 $j$ 求和（用 $\sum_j a_{ij}=1$）：

$$
\sum_{t=1}^{T-1}P(O,q_t=S_i\mid\lambda')+\eta_i=0
\ \Longrightarrow\
\eta_i=-\sum_{t=1}^{T-1}P(O,q_t=S_i\mid\lambda')
$$

代回：

$$
\hat a_{ij}=\frac{\sum_{t=1}^{T-1}P(O,q_t=S_i,q_{t+1}=S_j\mid\lambda')}{\sum_{t=1}^{T-1}P(O,q_t=S_i\mid\lambda')}
=\frac{\sum_{t=1}^{T-1}\xi_t(i,j)}{\sum_{t=1}^{T-1}\gamma_t(i)}
$$

$\square$

### 24.5.5 M 步（三）：$B$

第三项：

$$
\sum_{j=1}^{N}\sum_{t=1}^{T}\log b_j(o_t)\,P(O,q_t=S_j\mid\lambda')
$$

注意 $\log b_j(o_t)$ 只在 $o_t=v_k$ 时取 $\log b_j(k)$。引入示性函数 $\mathbb 1(o_t=v_k)$：

$$
\sum_{j=1}^{N}\sum_{t=1}^{T}\log b_j(o_t)P(O,q_t=S_j\mid\lambda')
=\sum_{j=1}^{N}\sum_{k=1}^{M}\Bigl(\sum_{t=1}^{T}\mathbb 1(o_t=v_k)P(O,q_t=S_j\mid\lambda')\Bigr)\log b_j(k)
$$

约束 $\sum_k b_j(k)=1$。对 $b_j(k)$ 求导置零，再对 $k$ 求和定乘子：

$$
\hat b_j(k)=\frac{\sum_{t=1}^{T}\mathbb 1(o_t=v_k)\,P(O,q_t=S_j\mid\lambda')}{\sum_{t=1}^{T}P(O,q_t=S_j\mid\lambda')}
=\frac{\sum_{t:\,o_t=v_k}\gamma_t(j)}{\sum_{t=1}^{T}\gamma_t(j)}
$$

$\square$

### 24.5.6 与 24.4 的对照

推导结果与 24.4 的「猜」完全一致。这不是巧合：**EM 的 M 步在指数族里就是「期望充分统计量」**，而 HMM 的三项都是多项分布，充分统计量就是计数。

## 24.6 一次迭代的数值（实跑）

用 22.5 的模型与 $O=(\text{R},\text{W},\text{R})$，跑一轮。

**$\gamma$ 表**

| $\gamma_t(i)$ | $i=1$ | $i=2$ | $i=3$ | 行和 |
|---|---|---|---|---|
| $t=1$ | 0.188223 | 0.322167 | 0.489610 | 1 |
| $t=2$ | 0.319311 | 0.415426 | 0.265263 | 1 |
| $t=3$ | 0.321538 | 0.272712 | 0.405750 | 1 |

**$\xi_t(i,j)$ 表**

$t=1$（$\xi_1$）：

| $i\backslash j$ | 1 | 2 | 3 |
|---|---|---|---|
| 1 | 0.103672 | 0.045155 | 0.039395 |
| 2 | 0.099525 | 0.180620 | 0.042022 |
| 3 | 0.116113 | 0.189651 | 0.183846 |

$t=2$（$\xi_2$）：

| $i\backslash j$ | 1 | 2 | 3 |
|---|---|---|---|
| 1 | 0.147829 | 0.047305 | 0.124176 |
| 2 | 0.127171 | 0.169562 | 0.118693 |
| 3 | 0.046537 | 0.055845 | 0.162881 |

两张表各元素之和都是 1（每个 $t$ 是一个概率分布）。

**重估**

$\hat\pi=\gamma_1=(0.1882,\ 0.3222,\ 0.4896)$

分母 $\sum_{t=1}^{2}\gamma_t=(0.507534,\ 0.737594,\ 0.754873)$，于是

$$
\hat A=\begin{bmatrix}
0.4955&0.1822&0.3223\\
0.3073&0.4748&0.2179\\
0.2155&0.3252&0.4593
\end{bmatrix}
$$

（例：$\hat a_{11}=(0.103672+0.147829)/0.507534=0.4955$）

$\hat B$ 的分子只累加 $o_t=R$ 的两个时刻（$t=1,3$）：

$$
\hat B=\begin{bmatrix}
0.6149&0.3851\\
0.5888&0.4112\\
0.7714&0.2286
\end{bmatrix}
$$

（例：$\hat b_1(R)=(0.188223+0.321538)/(0.188223+0.319311+0.321538)=0.6149$）

**似然变化**

$$
P(O\mid\lambda)=0.130218\ \longrightarrow\ P(O\mid\hat\lambda)=0.150463
$$

相对提升 $+15.6\%$。这一步是 EM 单调性的直接体现：**一次迭代后似然一定不降**。

> **补充** $T=3$、只有一条序列，这一轮之后模型已经强烈过拟合（$\hat b_3(R)=0.77$，$\hat\pi_3=0.49$）。这个例子的意义只是「演示方向与算术」，不是「训练结果」。真实训练要几十条序列、几十轮，且要留验证集看是否过拟合 —— 这一点与其他 ML 模型没有区别（[[06-2-监督学习·损失函数·梯度下降与数据集划分]]）。

## 24.7 收敛性与初值

**性质**（EM 的一般结论，HMM 继承）

1. 每轮 $P(O\mid\lambda)$ **单调不降**；
2. 收敛到**局部**极大（或鞍点），不保证全局；
3. 收敛速度在极值附近线性，通常比梯度法慢但每步更贵；
4. 参数会撞到**边界**（0 或 1）—— 因为 M 步是闭式解，没有正则项推它离开边界。

**初值的三个实务做法**

| 做法 | 说明 |
|---|---|
| 均匀初始化 | $\pi_i=1/N$，$a_{ij}=1/N$，$b_j(k)=1/M$。简单，但对称性容易让某些状态「永远分不到」 |
| 随机扰动 | 在均匀值上加小随机数，打破对称。最常用 |
| 用 Viterbi 的粗解初始化 | 对序列做 k-means 或按观测频率分段，得到初始状态标注，再按 24.3 数频数。收敛最快，但依赖初始分段质量 |

> **补充** 「状态永远分不到」在 HMM 里叫**状态死亡**：某个状态的 $\gamma_t(i)$ 全为 0，于是 $\hat\pi_i=0$、所有 $\hat a_{ij}$ 的分母为 0（0/0）。实现时必须加保护（分母加 $\epsilon$，或直接删掉该死状态）。这不是数值技巧，是算法的固有缺陷。

## 24.8 实现要点

**(1) 多序列**。真实训练是 $K$ 条独立序列 $\{O^{(k)}\}$，所有计数**跨序列求和**：

$$
\hat a_{ij}=\frac{\sum_{k}\sum_{t}\xi_t^{(k)}(i,j)}{\sum_{k}\sum_{t}\gamma_t^{(k)}(i)}
$$

**(2) 下溢**。$\alpha,\beta$ 都要在 log 域算，但 $\gamma=\alpha\beta/P$ 在 log 域要写成 $\exp(\log\alpha+\log\beta-\log P)$，且分母用 log-sum-exp。若用 Rabiner 的 scaling，则 $\alpha,\beta$ 各自带了不同的尺度因子，必须先还原尺度再相乘 —— 这是手写 BW 最常见的 bug。

**(3) 终止判据**。不要用「参数变化小于阈值」（参数可能来回摆），用**对数似然增量**：$\log P(O\mid\lambda^{(n+1)})-\log P(O\mid\lambda^{(n)})<\delta$（$\delta=10^{-4}$ 量级）。

**(4) 数值稳定性检查**。训练完必须验：每行 $A$、每行 $B$、$\pi$ 之和为 1（浮点误差 $10^{-12}$ 内）；没有负值；没有 NaN。出现过 NaN 基本就是「0/0」或「log(0)」。

## 订正

**1. 第 46 页（$A$ 的拉格朗日推导）写 $\eta=-P(O,q_t=S_i\mid\lambda')$，漏了对 $t$ 的求和。**

$\eta_i$ 是**常数**（对固定的 $i$），不可能带下标 $t$。正确是

$$
\eta_i=-\sum_{t=1}^{T-1}P(O,q_t=S_i\mid\lambda')
$$

同一页最后一行给出的 $\hat a_{ij}$ 公式是**对的**（分母确实带 $\sum_{t=1}^{T-1}$），所以只是中间一步漏写求和号。对照第 47 页（$B$ 那一项）写的是 $\eta=-\sum_{t=1}^{T}P(O,q_t=S_j\mid\lambda')$，带了求和 —— 两页不一致，以第 47 页为准。

**2. 第 39 页 $\gamma_t(i)=\sum_{j=1}^{N}\xi_t(i,j)$ 没限定 $t\le T-1$。**

$t=T$ 时 $\xi_T(i,j)$ 未定义（没有 $t+1$），此时只能用 $\gamma_T(i)=\alpha_T(i)/P(O\mid\lambda)$。若在代码里无条件套这个式子，最后一列 $\gamma$ 会全为 0，进而 $\hat b_j(k)$ 的分母漏掉最后一步。

**3. 第 41 页「初始状态 $\pi_i$ 的估计：为样本中初始状态为 $S_i$ 的频率（或 $t=1$ 时 $\gamma_1(i)$）」把两种情形的口径混在一行。**

前半句（数频数）只适用于**情形一**（已知状态序列）；后半句 $\gamma_1(i)$ 适用于**情形二**。无监督时没有「样本中初始状态」，$\pi$ 只能由 $\gamma_1$ 估计。这句话在本讲义里紧跟情形一的表格，很容易被读成「两种都行」。

**4. 第 40 页「状态 $S_i$ 转移的期望数 $\sum_{t=1}^{T-1}\gamma_t(i)$」的措辞有歧义。**

中文「$S_i$ 转移」既可读成「转移到 $S_i$」也可读成「从 $S_i$ 出发」。这里指后者（departures from $S_i$）。同一页另两行分别是「$S_i$ 出现的期望数 $\sum_{t=1}^{T}\gamma_t(i)$」与「$S_i$ 转移到 $S_j$ 的期望数 $\sum_{t=1}^{T-1}\xi_t(i,j)$」，三行只在上下限与符号上不同，抄错一个下标整套重估就崩。

---

系列：[[08-1-隐马尔可夫模型·定义与三个基本问题]] · [[08-2-概率计算·前向与后向算法]] · [[08-3-Viterbi解码·动态规划求最优路径]] · [[08-5-profile-HMM与HMMER]]
