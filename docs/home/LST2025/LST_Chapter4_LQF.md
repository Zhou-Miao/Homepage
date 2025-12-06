# Linear-quadratic forms

线性-二次型在空间计量经济学中非常重要，常用于推导估计量的渐近性质。设 $v_{n1}, \ldots, v_{nn}$ 是独立随机变量，具有有限的二阶矩。记 $V_n = (v_{n1}, \ldots, v_{nn})^\prime$，$A_n = (a_{ij})$ 是一个 $n \times n$ 非随机常数矩阵，$b_n = (b_{n1}, \ldots, b_{nn})^\prime$ 是一个 $n \times 1$ 常数列向量。那么：

$$
Q_n = V_n^\prime A_n V_n + b_n^\prime V_n
$$

被称为 $V_n$ 的**线性-二次型**。


## 2.1 基本矩公式

**引理 1（线性-二次型的矩）**：设 $A_n = [a_{ij}]$ 是 $n \times n$ 常数矩阵。假设 $v_{n1}, \ldots, v_{nn}$ 是 **i.i.d.** 的，满足：
- $\mathbb{E} v_i = 0$
- $\operatorname{Var}(v_i) = \sigma^2$
- 四阶矩 $\mu_4 = \mathbb{E}[v_i^4] < \infty$

则：

#### (1) $\mathbb{E}(V_n^\prime A_n V_n) = \sigma^2 \operatorname{tr}(A_n)$

*证明*：
$$
\begin{aligned}
\mathbb{E}(V_n^\prime A_n V_n) &= \mathbb{E}\left(\sum_{i=1}^n \sum_{j=1}^n a_{ij} v_i v_j\right) \\
&= \sum_{i=1}^n \sum_{j=1}^n a_{ij} \mathbb{E}[v_i v_j]
\end{aligned}
$$

由于 $v_i$ 独立且 $\mathbb{E} v_i = 0$：
- 当 $i = j$ 时：$\mathbb{E}[v_i v_j] = \mathbb{E}[v_i^2] = \sigma^2$
- 当 $i \neq j$ 时：$\mathbb{E}[v_i v_j] = \mathbb{E} v_i \cdot \mathbb{E} v_j = 0$

因此：
$$
\mathbb{E}(V_n^\prime A_n V_n) = \sum_{i=1}^n a_{ii} \sigma^2 = \sigma^2 \sum_{i=1}^n a_{ii} = \sigma^2 \operatorname{tr}(A_n)
$$

#### (2) $\mathbb{E}[(V_n^\prime A_n V_n)^2] = (\mu_4 - 3\sigma^4)\sum_{i=1}^n a_{ii}^2 + \sigma^4[\operatorname{tr}^2(A_n) + \operatorname{tr}(A_n A_n^\prime) + \operatorname{tr}(A_n^2)]$

*证明*：
<!-- 这是数学公式 -->
$$
\begin{aligned}
\mathbb{E}[(V_n^\prime A_n V_n)^2] &= \mathbb{E}\left[\left(\sum_{i=1}^n \sum_{j=1}^n a_{ij} v_i v_j\right)^2\right] \\
&= \mathbb{E}\left[\sum_{i=1}^n \sum_{j=1}^n \sum_{k=1}^n \sum_{l=1}^n a_{ij} a_{kl} v_i v_j v_k v_l\right]
\end{aligned}
$$
<!-- 公式结束 -->

由于 $v_i$ i.i.d. 且 $\mathbb{E} v_i = 0$，$\mathbb{E}[v_i v_j v_k v_l]$ 仅当索引配对时非零。有四种配对情况：

1. **情况 1**：$i = j = k = l$
   $$
   \mathbb{E}\left[\sum_{i=1}^n a_{ii}^2 v_i^4\right] = \mu_4 \sum_{i=1}^n a_{ii}^2
   $$

2. **情况 2**：$i = j \neq k = l$
   $$
   \mathbb{E}\left[\sum_{i=1}^n \sum_{k=1,k\neq i}^n a_{ii} a_{kk} v_i^2 v_k^2\right] = \sigma^4 \sum_{i=1}^n \sum_{k=1,k\neq i}^n a_{ii} a_{kk} = \sigma^4 (\operatorname{tr}^2(A_n) - \sum_{i=1}^n a_{ii}^2)
   $$

3. **情况 3**：$i = k \neq j = l$
   $$
   \mathbb{E}\left[\sum_{i=1}^n \sum_{j=1,j\neq i}^n a_{ij}^2 v_i^2 v_j^2\right] = \sigma^4 \sum_{i=1}^n \sum_{j=1,j\neq i}^n a_{ij}^2 = \sigma^4 (\operatorname{tr}(A_n A_n^\prime) - \sum_{i=1}^n a_{ii}^2)
   $$

4. **情况 4**：$i = l \neq j = k$
   $$
   \mathbb{E}\left[\sum_{i=1}^n \sum_{j=1,j\neq i}^n a_{ij} a_{ji} v_i^2 v_j^2\right] = \sigma^4 \sum_{i=1}^n \sum_{j=1,j\neq i}^n a_{ij} a_{ji} = \sigma^4 (\operatorname{tr}(A_n^2) - \sum_{i=1}^n a_{ii}^2)
   $$

将四种情况相加：
<!-- 这是数学公式 -->
$$
\begin{aligned}
\mathbb{E}[(V_n^\prime A_n V_n)^2] &= \mu_4 \sum_{i=1}^n a_{ii}^2 + \sigma^4 [(\operatorname{tr}^2(A_n) - \sum a_{ii}^2) \\
&\quad + (\operatorname{tr}(A_n A_n^\prime) - \sum a_{ii}^2) + (\operatorname{tr}(A_n^2) - \sum a_{ii}^2)] \\
&= (\mu_4 - 3\sigma^4) \sum_{i=1}^n a_{ii}^2 + \sigma^4 [\operatorname{tr}^2(A_n) + \operatorname{tr}(A_n A_n^\prime) + \operatorname{tr}(A_n^2)]
\end{aligned}
$$
<!-- 公式结束 -->

#### (3) $\operatorname{Var}(V_n^\prime A_n V_n) = (\mu_4 - 3\sigma^4)\sum_{i=1}^n a_{ii}^2 + \sigma^4[\operatorname{tr}(A_n A_n^\prime) + \operatorname{tr}(A_n^2)]$

*证明*：由方差定义 $\operatorname{Var}(X) = \mathbb{E}[X^2] - (\mathbb{E}X)^2$，结合 (1) 和 (2)：
<!-- 这是数学公式 -->
$$
\begin{aligned}
\operatorname{Var}(V_n^\prime A_n V_n) &= \mathbb{E}[(V_n^\prime A_n V_n)^2] - [\mathbb{E}(V_n^\prime A_n V_n)]^2 \\
& = (\mu_4 - 3\sigma^4)\sum_{i=1}^n a_{ii}^2 + \sigma^4[\operatorname{tr}^2(A_n) + \operatorname{tr}(A_n A_n^\prime) + \operatorname{tr}(A_n^2)] - \sigma^4 \operatorname{tr}^2(A_n) \\
& = (\mu_4 - 3\sigma^4)\sum_{i=1}^n a_{ii}^2 + \sigma^4[\operatorname{tr}(A_n A_n^\prime) + \operatorname{tr}(A_n^2)]
\end{aligned}
$$
<!-- 公式结束 -->

**特殊情况**：

- 当 $v_i \sim N(0, \sigma^2)$ 时，$\mu_4 = 3\sigma^4$，方差简化为 $\sigma^4[\operatorname{tr}(A_n A_n^\prime) + \operatorname{tr}(A_n^2)]$。

- 当 $A_n$ 对称且幂等时，$\operatorname{tr}(A_n^2) = \operatorname{tr}(A_n)$，$\operatorname{tr}(A_n A_n^\prime) = \operatorname{tr}(A_n^2)$。



## 2.2 线性-二次型的大数定律

**引理 2（线性-二次型的弱大数定律）**：假设 $\{A_n\}$ 的行和与列和一致有界，$v_{n1}, \ldots, v_{nn}$ 是 i.i.d. 的，满足 $\mathbb{E} v_i = 0$ 且 $\operatorname{Var}(v_i) = \sigma^2 < \infty$。那么：

1. $\mathbb{E}(V_n^\prime A_n V_n) = O(n)$
2. $\operatorname{Var}(V_n^\prime A_n V_n) = O(n)$
3. $V_n^\prime A_n V_n = O_p(n)$
4. $\frac{1}{n} V_n^\prime A_n V_n - \frac{1}{n} \mathbb{E}(V_n^\prime A_n V_n) = O_p(n^{-1/2})$

*完整证明*：

(1) 由引理 1 (1)：
$$
\mathbb{E}(V_n^\prime A_n V_n) = \sigma^2 \operatorname{tr}(A_n)
$$
由于 $\{A_n\}$ 的行和一致有界，设 $\|A_n\|_\infty \leq M$，则 $|a_{ii}| \leq M$，所以：
$$
|\operatorname{tr}(A_n)| = \left|\sum_{i=1}^n a_{ii}\right| \leq \sum_{i=1}^n |a_{ii}| \leq nM
$$
因此 $\mathbb{E}(V_n^\prime A_n V_n) = O(n)$。

(2) 由引理 1 (3)：
$$
\operatorname{Var}(V_n^\prime A_n V_n) = (\mu_4 - 3\sigma^4)\sum_{i=1}^n a_{ii}^2 + \sigma^4[\operatorname{tr}(A_n A_n^\prime) + \operatorname{tr}(A_n^2)]
$$

对于第一项：$\sum_{i=1}^n a_{ii}^2 \leq n \max_i a_{ii}^2 \leq nM^2 = O(n)$

对于第二项：$\operatorname{tr}(A_n A_n^\prime) = \sum_{i=1}^n \sum_{j=1}^n a_{ij}^2$

由于 $\{A_n\}$ 的行和一致有界，设 $\sum_{j=1}^n |a_{ij}| \leq M$，由 Cauchy-Schwarz 不等式：
$$
\sum_{j=1}^n a_{ij}^2 \geq \frac{(\sum_{j=1}^n |a_{ij}|)^2}{n} \geq 0
$$
且
$$
\sum_{j=1}^n a_{ij}^2 \leq (\sum_{j=1}^n |a_{ij}|)^2 \leq M^2
$$
因此：
$$
\operatorname{tr}(A_n A_n^\prime) = \sum_{i=1}^n \sum_{j=1}^n a_{ij}^2 \leq \sum_{i=1}^n M^2 = nM^2 = O(n)
$$

同理，$\operatorname{tr}(A_n^2) = O(n)$。

综上，$\operatorname{Var}(V_n^\prime A_n V_n) = O(n)$。

(3) 由 Chebyshev 不等式：
$$
P\left(|V_n^\prime A_n V_n - \mathbb{E}(V_n^\prime A_n V_n)| > \epsilon n\right) \leq \frac{\operatorname{Var}(V_n^\prime A_n V_n)}{(\epsilon n)^2} = \frac{O(n)}{n^2} \to 0
$$
所以 $V_n^\prime A_n V_n = \mathbb{E}(V_n^\prime A_n V_n) + O_p(n) = O_p(n)$。

(4) 由 $\operatorname{Var}(V_n^\prime A_n V_n) = O(n)$，有：
$$
\operatorname{Var}\left(\frac{1}{n} V_n^\prime A_n V_n\right) = \frac{1}{n^2} \operatorname{Var}(V_n^\prime A_n V_n) = O\left(\frac{1}{n}\right)
$$
因此：
$$
\frac{1}{n} V_n^\prime A_n V_n - \frac{1}{n} \mathbb{E}(V_n^\prime A_n V_n) = O_p(n^{-1/2})
$$

---

## 2.3 鞅差分解

线性-二次型的中心极限定理依赖于**鞅差序列**的中心极限定理。

**定义 3（适应随机序列）**：设 $\{z_{nt}: t = 1, \ldots, n\}$ 是随机变量的三角阵列，$\{\mathcal{J}_{nt}\}$ 是 $\sigma$-域的三角阵列，满足 $\mathcal{J}_{n,t-1} \subseteq \mathcal{J}_{nt}$ 对所有 $n$ 和 $t$。如果 $z_{nt}$ 关于 $\mathcal{J}_{nt}$ 可测，则 $\{z_{nt}, \mathcal{J}_{nt}: t = 1, \ldots, n\}$ 称为**适应随机序列 (adapted stochastic sequence)**。

- $\sigma$-域 $\mathcal{J}_{nt}$ 可以看作由 $z_{nt}$ 和一些其他随机变量的当前和过去历史生成的 $\sigma$-代数。

**定义 4（鞅差阵列）**：设 $\{y_{nt}, \mathcal{J}_{nt}\}$ 是适应随机序列。那么 $\{y_{nt}, \mathcal{J}_{nt}\}$ 是**鞅差阵列（MDA）** 当且仅当 $\mathbb{E}(y_{nt} | \mathcal{J}_{n,t-1}) = 0$ 对所有 $n$ 和 $t \geq 2$ 成立。

**引理 5（线性-二次型的鞅差分解）**：假设 $v_{ni}$ 是 i.i.d. 的，满足 $\mathbb{E} v_{ni} = 0$ 和有限方差。对于线性-二次型 $Q_n = V_n^\prime A_n V_n + b_n^\prime V_n$，那么 $ Q_n - \mathbb{E} Q_n$ 可以写成一个鞅差阵列的和。

*证明*：
将 $Q_n - \mathbb{E} Q_n$ 重写为：
<!-- 这是数学公式 -->
$$
\begin{aligned}
Q_n - \mathbb{E} Q_n &= \sum_{i=1}^n b_i v_i + \sum_{i=1}^n \sum_{j=1}^n a_{ij} v_i v_j - \sigma^2 \operatorname{tr}(A_n) \\
&= \sum_{i=1}^n b_i v_i + \sum_{i=1}^n a_{ii} (v_i^2 - \sigma^2) + 2\sum_{i=1}^n \sum_{j=1}^{i-1} a_{ij} v_i v_j
\end{aligned}
$$
<!-- 公式结束 -->

现在，将上式重新组织为求和形式。定义：
$$
Z_{ni} = b_i v_i + a_{ii}(v_i^2 - \sigma^2) + 2v_i \sum_{j=1}^{i-1} a_{ij} v_j
$$
则：
$$
Q_n - \mathbb{E} Q_n = \sum_{i=1}^n Z_{ni}
$$

定义 $\sigma$-域 $\mathcal{J}_{ni} = \sigma\{v_1, \ldots, v_i\}$（由 $v_1, \ldots, v_i$ 生成的 $\sigma$-代数）。由于 $v_i$ 是 i.i.d. 的且均值为零，有：
<!-- 这是数学公式 -->
$$
\begin{aligned}
\mathbb{E}(Z_{ni} | \mathcal{J}_{n,i-1}) &= \mathbb{E}(b_i v_i | \mathcal{J}_{n,i-1}) + \mathbb{E}(a_{ii}(v_i^2 - \sigma^2) | \mathcal{J}_{n,i-1}) + \mathbb{E}(2v_i \sum_{j=1}^{i-1} a_{ij} v_j | \mathcal{J}_{n,i-1}) \\
&= b_i \mathbb{E}(v_i) + a_{ii}(\mathbb{E}(v_i^2) - \sigma^2) + 2\sum_{j=1}^{i-1} a_{ij} v_j \mathbb{E}(v_i) \\
&= 0 + a_{ii}(\sigma^2 - \sigma^2) + 0 \\
&= 0
\end{aligned}
$$
<!-- 公式结束 -->

因此 $\{(Z_{ni}, \mathcal{J}_{ni})\}$ 构成一个鞅差三角阵列。

---

## 2.4 中心极限定理

**定理 6（鞅差中心极限定理）**：设 $\{\xi_{Tt}, \mathcal{J}_{T,t}\}$ 是鞅差阵列。如果当 $T \to \infty$ 时：

(i) $\sum_{t=1}^T v_{Tt} \overset{p}{\to} \sigma^2$，其中 $\sigma^2 > 0$，$v_{Tt} = \mathbb{E}(\xi_{Tt}^2 | \mathcal{J}_{T,t-1})$ 是条件方差；

(ii) 对任意 $\epsilon > 0$，$\sum_{t=1}^T \mathbb{E}[\xi_{T_t}^2 I(|\xi_{T_t}| > \epsilon) | \mathcal{J}_{T,t-1}] \overset{p}{\to} 0$（Lindeberg 条件）；

则：
$$
\xi_{T1} + \cdots + \xi_{TT} \overset{d}{\to} N(0, \sigma^2).
$$

**定理 7（Kelejian and Prucha, 2001）**：设 $V_n = (v_{1,n}, \ldots, v_{n,n})$ 的元素是独立的，满足 $\mathbb{E} v_{i,n} = 0$。$A_n = (a_{ij,n})$ 是非随机对称矩阵，其列和一致有界。$b_n = (b_{1,n}, \ldots, b_{n,n})$ 是非随机向量，满足 $\sup_n \frac{1}{n} \sum_{i=1}^n |b_{i,n}|^{2+\eta_1} < \infty$ 对某个 $\eta_1 > 0$。假设 $\sup_{i,n} \mathbb{E}|v_{i,n}|^{4+\eta_2} < \infty$ 对某个 $\eta_2 > 0$ 成立。记 $Q_n \equiv V_n^\prime A_n V_n + b_n^\prime V_n$，$\sigma_{Q_n} \equiv [\operatorname{Var}(Q_n)]^{1/2}$。如果 $n^{-1} \sigma_{Q_n}^2 \geq c$ 对某个 $c > 0$ 对所有 $n$ 成立，那么：
$$
\frac{Q_n - \mathbb{E} Q_n}{\operatorname{std}(Q_n)} \overset{d}{\to} N(0, 1)
$$

*证明思路*：
1. 首先通过引理 5 将 $Q_n$ 分解为鞅差和
2. 验证定理 6 的条件 (i) 和 (ii)
3. 应用鞅差中心极限定理得到渐近正态性

条件 $n^{-1} \sigma_{Q_n}^2 \geq c$ 确保方差以 $n$ 的速度增长，避免了退化的极限分布。


