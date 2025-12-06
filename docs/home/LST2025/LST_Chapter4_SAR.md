# 空间计量经济学 （Spatial Econometrics）


## 1. Introduction
由于空间溢出效应、网络效应与共同冲击，经济金融、社会科学与自然科学（如气象、海洋、医学）中的数据常常存在**截面相关性**。传统上对截面数据常假设独立甚至同分布，因此自然的问题是：如何建模截面相关性？在计量与统计中，主要有以下几种方法：

- **空间计量模型**：假设个体在空间或网络中存在交互或溢出效应，适用于建模*弱空间相关*，即相关性主要体现在*邻近个体之间*。
- **因子模型**：主要用于建模共同冲击（强相关），如金融危机期间多数股票价格同时下跌。
- **协方差结构**：直接建模协方差矩阵，如股票收益的协方差结构。

### 1.1 常见的空间计量模型
以下是几种常见的空间计量模型：

| 模型名称 | 模型形式 | 说明 |
|---|---|---|
| **纯SAR模型** | $Y_n = \lambda W_n Y_n + V_n$ | 仅包含空间滞后项，无外生变量 |
| **SAR模型** | $Y_n = \lambda W_n Y_n + X_n \beta + V_n$ | 包含外生变量 $X_n$，最常用 |
| **高阶SAR模型** | $Y_n = \sum_{l=1}^p \lambda_l W_{ln} Y_n + X_n \beta + V_n$ | 使用多个空间权重矩阵 |
| **空间误差模型（SE）** | $Y_n = X_n \beta + U_n,\quad U_n = \lambda W_n U_n + V_n$ | 误差项具有空间自相关 |
| **SARAR模型** | $Y_n = \lambda W_n Y_n + U_n,\quad U_n = \rho M_n U_n + V_n$ | 同时包含空间滞后与空间误差 |
| **带Durbin项的SAR模型** | $Y_n = \lambda W_n Y_n + l_n \alpha + X_n \beta_1 + W_n X_n \beta_2 + V_n$ | 引入 $W_n X_n$ 捕捉外生邻里效应 |
| **带外生邻居效应的回归模型** | $Y_n = X_n \beta + W_n X_n \gamma + V_n$ | 仅通过 $W_n X_n$ 引入空间交互 |


### 1.2 空间权重矩阵

空间权重矩阵 $W_n$ 是一个 $n\times n$ 矩阵，用于量化空间或网络中单位之间的连接强度或距离关系。其一般性质与构造方法如下：

#### 1.2.1 基本性质
- 对角线元素通常为 0：$w_{ii} = 0$，表示单位对自身无空间影响。
- 非对角线元素 $w_{ij}$ 表示单位 $j$ 对单位 $i$ 的影响强度。
- 矩阵通常**非负**：$w_{ij} \geq 0$。

#### 1.2.2 常见构造方法
#### 1. 基于距离的函数
设 $d_{ij}$ 为两单位间的距离，定义：

$$
h(d_{ij}) = d_{ij}^{-\alpha} \quad \text{或} \quad h(d_{ij}) = \exp(-\alpha d_{ij}), \quad \alpha > 0
$$

则权重为：

$$
w_{ij} = \frac{h(d_{ij})}{\sum_{k=1}^n h(d_{ik})}
$$

此时 $W$ 经过行标准化（row-normalization）。

#### 2. 邻接矩阵（0-1矩阵）
- 若 $i$ 与 $j$ 相邻（如地理接壤、网络连接），则 $w_{ij} = 1$，否则为 0。
- 可进一步行标准化：若单位 $i$ 有 $n_i$ 个邻居，则：
$$
w_{ij} = \frac{1}{n_i} \quad (\text{若 } j \text{ 是 } i \text{ 的邻居})
$$

#### 3. 基于群组结构
假设样本分为 $R$ 组，每组 $m$ 个单位，组内权重均匀分配：

$$
B_m = \frac{1}{m-1}(l_m l_m' - I_m), \quad W_n = I_R \otimes B_m
$$

其中 $l_m$ 是全1向量，$\otimes$ 为Kronecker积。

#### 1.2.3 行标准化的意义与问题
- **意义**：使得 $\sum_j w_{ij} = 1$，此时 $W_n Y_n$ 可解释为邻居观测值的加权平均。
- **问题**：若 $W_n$ 行标准化且 $X_n$ 包含截距项，则 $W_n X_n$ 与 $X_n$ 可能产生多重共线性。此时应去除 $W_n X_n$ 中的截距列。

#### 1.2.4 注意事项
- $W_n$ 不一定对称，尤其在行标准化后。
- 为保证估计理论的大样本性质，常假设 $\{W_n\}$ 的行和与列和一致有界（即不存在"主导单位"）。
- 若存在主导单位（如国际贸易中的大国），需采用 [Lee, Yang & Yu (2022) 等：QML and Efficient GMM Estimation of Spatial Autoregressive Models with Dominant (Popular) Units](https://www.tandfonline.com/doi/full/10.1080/07350015.2022.2041424#abstract) 扩展方法。


### 1.3 空间计量中的矩阵范数应用

在空间计量中，矩阵范数常用于分析空间权重矩阵 $W_n$ 的性质。一个重要结果是：

当 $||\lambda W_n||_\infty < 1$ 时（例如 $||W_n||_\infty = 1$ 且 $|\lambda| < 1$），有：
$$
(I_n - \lambda W_n)^{-1} = \sum_{j=0}^\infty (\lambda W_n)^j = I_n + \lambda W_n + \lambda^2 W_n^2 + \cdots
$$

这个展开在解释 SAR 模型中的溢出效应时非常重要：
- $\lambda W_n$ 表示直接邻居的影响
- $\lambda^2 W_n^2$ 表示二阶邻居的影响（邻居的邻居）
- 以此类推，高阶项表示更远距离的间接影响

对于 SAR 模型 $Y_n = \lambda W_n Y_n + X_n \beta + V_n$，其简化形式为：
$$
Y_n = \sum_{j=0}^\infty \lambda^j W_n^j (X_n \beta + V_n)
$$
显示了外生变量 $X_n$ 和误差 $V_n$ 通过空间结构的传播机制。

此外，我们有级数截断的误差界：

$$
\left\|\sum_{j=m}^{\infty} \lambda^j W_n^j\right\|_{\infty} \leq \sum_{j=m}^{\infty} \|\lambda W_n^j\|_{\infty} \leq \sum_{j=m}^{\infty} \|\lambda W_n\|^j_{\infty}\leq \sum_{j=m}^{\infty} |\lambda|^j = \frac{|\lambda|^m}{1-|\lambda|}
$$

当 $m$ 足够大时，$\frac{|\lambda|^m}{1-|\lambda|}$ 很小。在实际应用中，常使用有限项近似：

$$
(I_n - \lambda W_n)^{-1} \approx I_n + \lambda W_n + \lambda^2 W_n^2 + \cdots + \lambda^m W_n^m
$$

其中 $m$ 的选择使得 $\left\|\sum_{j=m+1}^\infty \lambda^j W_n^j\right\|$ 足够小。


#### 1.3.1 空间权重矩阵的规范性条件
在空间计量理论中，常假设空间权重矩阵 $\{W_n\}$ 和 $\{S_n^{-1}\}$（其中 $S_n = I_n - \lambda W_n$）在行和与列和上一致有界。等价地说，假设 $\{\|W_n\|_1\}$ 和 $\{\|W_n\|_\infty\}$ 是有界序列，同样 $S_n^{-1}$ 的最大行和与列和范数也有界。

这一假设的一个重要性质是：如果 $\{A_n\}$ 和 $\{B_n\}$ 在行和（列和）上一致有界，那么 $\{A_n B_n\}$ 也在行和（列和）上一致有界。这源于矩阵范数的次可乘性，例如：当 $\|A_n\|_\infty \leq c$ 且 $\|B_n\|_\infty \leq c$ 时，
<!-- 这是数学公式 -->
$$
\|A_n B_n\|_\infty \leq \|A_n\|_\infty \|B_n\|_\infty \leq c^2.
$$
<!-- 公式结束 -->

**经济含义**：这一假设意味着不存在"主导"的空间单位，即所有空间单位的影响都是有界的。这在许多现实经济应用中是一个合理的假设，但在某些情况下（如国际贸易中美国和中国这样的大国）可能不成立。当存在主导空间单位时，需要参考 [Lee, Yang & Yu (2022) 等：QML and Efficient GMM Estimation of Spatial Autoregressive Models with Dominant (Popular) Units](https://www.tandfonline.com/doi/full/10.1080/07350015.2022.2041424#abstract) 等扩展方法。


---

## 2. SAR 模型的估计

假设 SAR 模型中空间权重矩阵 $W_n$ 与 $X_n$ 是非随机的，

$$
Y_n = \lambda_0 W_n Y_n + X_n \beta_0 + V_n.
$$


记 $S_n(\lambda) \equiv I_n - \lambda W_n, S_n \equiv I_n - \lambda_0 W_n$。 假设 $S_n(\lambda)$ 对于参数空间中的所有 $\lambda$ 都是可逆的，则
$$
Y_n = (I_n - \lambda_0 W_n)^{-1} ( X_n \beta_0 + V_n) = S_n(\lambda)^{-1} ( X_n \beta_0 + V_n).
$$
因此，$W_n Y_n = W_n S_n(\lambda)^{-1} ( X_n \beta_0 + V_n) \equiv G_n ( X_n \beta_0 + V_n)$. 从而我们有
$$
\mathbb{E}[(W_n Y_n)'V_n] = \mathbb{E} [( X_n \beta_0 + V_n)'G_n'V_n] = \mathbb{E} [V_n' G_n'V_n ] = \sigma^2 \operatorname{tr}(G_n).
$$
也就是说，$W_n Y_n$ 是一个内生变量 (endogenous variable)。


当存在内生变量时，通常有三种方法来处理：

- 工具变量/两阶段最小二乘；
- 联立方程估计： 如果内生性可以被建模，可以同时估计方程组；
- 控制函数 (control function)。

### 2.1 两阶段最小二乘 (2SLS)

首先我们介绍 2SLS 估计。

???+ question "一般地，工具变量应该满足什么条件？"

    - **相关性**：工具变量与内生变量相关；
    - **外生性**：工具变量与误差项不相关。

1. 假设 $\| \lambda_0 W_n\|_\infty < 1$，则 $S_n^{-1} \equiv (I_n - \lambda_0 W_n)^{-1} = \sum_{l=0}^\infty \lambda_0^l W_n^l$，并且
$$
W_n Y_n = W_n S_n^{-1} (X_n \beta_0 + V_n) = \sum_{l=0}^\infty \lambda_0^l W_n^{l+1} (X_n \beta_0 + V_n).
$$
因此，$W_n X_n, W_n X_n^2,\ldots$ 通常与 $W_n Y_n$ 相关。注意，这里的相关性是从样本的角度来说的，即当 $n\rightarrow \infty$时， $\frac{1}{n} (W_n X_n)'(W_n Y_n)$ 不会衰减。

2. 因为 $W_n$ 与 $X_n$ 是非随机的，$\mathbb{E}(V_n)=0$，所以 $W_n X_n, W_n X_n^2,\ldots$ 与 $V_n$ 不相关，即 $\mathbb{E}(X_n'V_n)=0$。

因此，$W_nX_n,W_n^2X_n$ 可以作为 IV。通常，我们使用 $[W_nX_n ~ ~ X_n]$ 或者 $[W_nX_n ~ ~ W_n^2X_n  ~ ~ X_n]$ 作为工具变量矩阵，记作 $Q_n= (q_{1n},\ldots,q_{nn})'$。我们也可以假设 $X_n$ 是随机的，但需满足 $\mathbb{E}(V_n|X_n) = 0$。

#### 2.1.1 两阶段最小二乘估计量
记 $Z_n = (W_n Y_n, X_n)$，$\theta = (\lambda, \beta')'$。SAR 模型可重写为 $Y_n = Z_n \theta + V_n$。使用工具变量 $Q_n$ 的 2SLS 估计量为：

$$
\hat{\theta}_{2sls} = [Z_n' Q_n (Q_n' Q_n)^{-1} Q_n' Z_n]^{-1} Z_n' Q_n (Q_n' Q_n)^{-1} Q_n' Y_n.
\tag{1}
$$

记投影矩阵 $P_Q = Q_n (Q_n' Q_n)^{-1} Q_n'$，则

$$
\hat{\theta}_{2sls} = [Z_n' P_Q Z_n]^{-1} Z_n' P_Q Y_n.
$$

**注意**：当 $W_n$ 行标准化时，$W_n \iota_n = \iota_n$。因此，若 $X_n$ 的第一列为 $\iota_n$，则应将其从 $W_n X_n$ 中剔除（设 $X_n$ 的第一列是 $\iota_n$）。记 $G_n \equiv W_n S_n^{-1}$。

$$
\mathbb{E} Z_n' Q_n = [W_n S_n^{-1} X_n \beta_0 ~ ~ X_n]' Q \equiv [G_n X_n \beta_0 ~ ~ X_n]' Q.
$$



#### 2.1.2 渐近性质的假设

**假设 1**：$x_{i,n}$ 的元素一致有界；$W_n, S_n^{-1}$ 的行范数与列范数一致有界。

**假设 2**：$v_i$ 是 i.i.d. $(0, \sigma^2)$。

**假设 3**：
(i) $M_{ZQ} \equiv \lim_{n \to \infty} \frac{1}{n} (G_n X_n \beta_0, X_n)' Q_n$ 存在且满秩；
(ii) $\Sigma_{QQ} \equiv \lim_{n \to \infty} \frac{1}{n} Q_n' Q_n$ 存在且正定；
(iii) $\sup_n \left( \frac{1}{n} \sum_{i=1}^n \| q_{in} \|^{2+\delta} \right) < \infty$ 对某个 $\delta > 0$。


#### 2.1.3 2SLS 的渐近正态性

**定理 1 （2SLS 的渐近正态性）**：在假设 1–3 下，

$$
  \sqrt{n} (\hat{\theta}_{2sls} - \theta_0) \overset{d}{\to} N\left( 0, \sigma^2 (M_{ZQ} \Sigma_{QQ}^{-1} M_{ZQ}')^{-1} \right).
$$

*证明*：由于

$$
\sqrt{n} (\hat{\theta}_{2sls} - \theta_0) = \left[ \frac{1}{n} Z_n' Q_n \left( \frac{1}{n} Q_n' Q_n \right)^{-1} \frac{1}{n} Q_n' Z_n \right]^{-1} \frac{1}{n} Z_n' Q_n \left( \frac{1}{n} Q_n' Q_n \right)^{-1} \frac{1}{\sqrt{n}} Q_n' V_n,
$$

在假设 3 下，由 Lyapunov CLT 与 Cramer–Wold 方法可得

$$
\frac{1}{\sqrt{n}} Q_n' V_n = \frac{1}{\sqrt{n}} \sum_{i=1}^n q_{i,n} v_{i,n} \overset{d}{\to} N(0, \sigma^2 \Sigma_{QQ}).
$$

若 $\frac{1}{n} Z_n' Q_n \overset{p}{\to} M_{ZQ}$，则根据 Slutsky 引理，

$$
\sqrt{n} (\hat{\theta}_{2sls} - \theta_0) \overset{d}{\to} (M_{ZQ} \Sigma_{QQ}^{-1} M_{ZQ}')^{-1} M_{ZQ} \Sigma_{QQ}^{-1} \cdot N(0, \sigma^2 \Sigma_{QQ}) = N(0, \sigma^2 (M_{ZQ} \Sigma_{QQ}^{-1} M_{ZQ}')^{-1}).
$$

因此只需证明 $\frac{1}{n} Z_n' Q_n \overset{p}{\to} M_{ZQ}$。

由于 $Z_n' Q_n$ 的第二部分 $X_n' Q$ 是非随机的，故只需证明

$$
\frac{1}{n} (W_n Y_n - G_n X_n \beta_0)' Q_n = \frac{1}{n} (G_n V_n)' Q_n = o_p(1).
$$

因 $\mathbb{E} \frac{1}{n} (G_n V_n)' Q_n = 0$，只需证

$$
\operatorname{Var} \left[ \frac{1}{n} (G_n V_n)' Q_n \right] = \frac{\sigma^2}{n^2} Q_n' G_n G_n' Q_n = o(1).
$$

由半正定矩阵 $\Sigma$ 的性质，其所有主子式非负，即 $\Sigma_{ii} \Sigma_{jj} - \Sigma_{ij}^2 \ge 0$，故只需证 $\frac{1}{n^2} Q_n' G_n G_n' Q_n$ 的所有对角元为 $o(1)$：

$$
\frac{1}{n^2} (Q_{.i,n}' G_n G_n' Q_{.i,n}) = o(1),
$$

其中 $Q_{.i,n}$ 为 $Q_n$ 的第 $i$ 列。

由于 $\sup_n (\| G_n \|_1 + \| G_n \|_\infty) < \infty$，有

$$
\sup_n \max \operatorname{eig}(G_n G_n') \le \sup_n \| G_n G_n' \|_\infty \le \sup_n \| G_n \|_1 \cdot \| G_n \|_\infty < \infty.
$$

又 $\lim_{n \to \infty} \frac{1}{n} Q_n' Q_n$ 存在，故

$$
0 \le \frac{1}{n^2} (Q_{.i,n}' G_n G_n' Q_{.i,n}) \le \max \operatorname{eig}(G_n G_n') \cdot \frac{1}{n^2} (Q_{.i,n}' Q_{.i,n}) = O(1) \cdot O(n^{-1}) = O_p(n^{-1}).
$$

因此 $\frac{1}{n} Z_n' Q_n \overset{p}{\to} M_{ZQ}$，定理得证。


#### 2.1.4 最优工具变量 (Best IV)

[Lee (2003): Best Spatial Two-Stage Least Squares Estimators for a Spatial Autoregressive Model with Autoregressive Disturbances.](https://www.tandfonline.com/doi/full/10.1081/ETC-120025891) 提出了 SAR 模型的最优 IV 估计。因为 $\mathbb{E} W_n Y_n = W_n S_n^{-1} X_n \beta_0 = G_n X_n \beta_0$，所以 $W_n Y_n$ 的最优工具变量应为 $G_n X_n \beta_0$。

严格论证如下：回顾 $\mathbb{E} Z_n' Q_n \equiv [G_n X_n \beta_0 ~ ~ X_n] Q$。由定理 1，

$$
\operatorname{avar}[\sqrt{n} (\hat{\theta}_{2sls} - \theta_0)] = \sigma^2 (M_{ZQ} \Sigma_{QQ}^{-1} M_{ZQ}')^{-1} = \lim_{n \to \infty} \sigma^2 \left\{ \frac{1}{n} \mathbb{E} Z_n' Q_n \left( \frac{1}{n} Q_n' Q_n \right)^{-1} \frac{1}{n} Q_n' \mathbb{E} Z_n \right\}^{-1}.
$$

当 $Q_n = \mathbb{E} Z_n$ 时，

$$
\operatorname{avar}[\sqrt{n} (\hat{\theta}_{2sls} - \theta_0)] = \lim_{n \to \infty} \sigma^2 \left\{ \frac{1}{n} \mathbb{E} Z_n' \mathbb{E} Z_n \right\}^{-1}.
$$

因为 $P_Q \equiv Q_n (Q_n' Q_n)^{-1} Q_n'$ 是投影矩阵，有

$$
\mathbb{E} Z_n' Q_n (Q_n' Q_n)^{-1} Q_n' \mathbb{E} Z_n \le \mathbb{E} Z_n' \mathbb{E} Z_n,
$$

其中 $A \le B$ 表示 $B - A$ 半正定。因此，

$$
\lim_{n \to \infty} n \sigma^2 \left\{ \mathbb{E} Z_n' Q_n (Q_n' Q_n)^{-1} Q_n' \mathbb{E} Z_n \right\}^{-1} \ge \lim_{n \to \infty} n \sigma^2 \left\{ \mathbb{E} Z_n' \mathbb{E} Z_n \right\}^{-1}.
$$

换言之，在 $v_i$ i.i.d. $(0, \sigma^2)$ 时，**$\mathbb{E} Z_n$ 是最优工具变量**，且 **$G_n X_n \beta_0$ 是 $W_n Y_n$ 的最优 IV**。

#### 2.1.5 可行的最优 IV 估计

实际中 $\beta_0$ 未知，可行的最优 IV 估计可按以下两步进行：

1. 以 $W_n X_n$ 和 $W_n^2 X_n$ 作为 IV，用 (1) 得到 2SLS 估计量 $\hat{\beta}_{2sls}$；
2. 以 $\hat{G}_n X_n \hat{\beta}_{2sls}$ 作为 $W_n Y_n$ 的 IV，得到可行最优 IV 估计量：

$$
\hat{\theta}_{best} = (\hat{Z}_n' Q_n)^{-1} \hat{Z}_n' Y_n,
$$

其中 $\hat{Z}_n \equiv [\hat{G}_n X_n \hat{\beta}_{2sls} ~ ~ X_n]$，$\hat{G}_n = W_n (I_n - \hat{\lambda}_{2sls} W_n)^{-1}$。

#### 2.1.6 对应代码实现

``` R title="2SLS for SAR model"
####################################################
# Three main functions for SAR model 2SLS estimation
####################################################

# IV: (X, WX*)
SAR_2SLS1 <- function(Y, X, W){
  n = length(Y)
  
  X_star <- X[, -1, drop=FALSE]
  
  # IV matrix
  WX_star <- W %*% X_star
  Q <- cbind(X, WX_star)
  
  # endogenous variable
  WY <- W %*% Y
  Z <- cbind(WY, X)
  
  # First stage
  P_Q <- Q %*% solve(t(Q) %*% Q) %*% t(Q)
  Z_hat <- P_Q %*% Z
  
  # Second stage
  theta_hat <- solve(t(Z_hat) %*% Z_hat) %*% t(Z_hat) %*% Y
  
  return(list(
    lambda = theta_hat[1],
    beta = theta_hat[-1]
  ))
}


# IV: (X, WX*, W^2X*)
SAR_2SLS2 <- function(Y, X, W){
  n = length(Y)
  
  X_star <- X[, -1, drop=FALSE]
  
  # IV matrix
  WX_star <- W %*% X_star
  W2X_star <- W %*% WX_star
  Q <- cbind(X, WX_star, W2X_star)
  
  # endogenous variable
  WY <- W %*% Y
  Z <- cbind(WY, X)
  
  # First stage
  P_Q <- Q %*% solve(t(Q) %*% Q) %*% t(Q)
  Z_hat <- P_Q %*% Z
  
  # Second stage
  theta_hat <- solve(t(Z_hat) %*% Z_hat) %*% t(Z_hat) %*% Y
  
  return(list(
    lambda = theta_hat[1],
    beta = theta_hat[-1]
  ))
}



# feasible best IV
SAR_BIV <- function(Y, X, W){
  n <- length(Y)
  
  X_star <- X[, -1, drop=FALSE]
  
  # endogenous variable
  WY <- W %*% Y
  Z <- cbind(WY, X)
  
  # step 1: initial estimate via SAR_2SLS2
  init_fit <- SAR_2SLS2(Y, X, W)
  lambda_init <- init_fit$lambda
  beta_init <- init_fit$beta
  
  # step 2: construct best IV
  S_inv <- solve(diag(n) - lambda_init * W)
  G_n <- W %*% S_inv
  GX_beta <- G_n %*% X %*% beta_init
  Q <- cbind(GX_beta, X) # best IV matrix
  
  # # best IV estimation
  # # first stage
  # P_Q <- Q %*% solve(t(Q) %*% Q) %*% t(Q)
  # Z_hat <- P_Q %*% Z
  # # second stage
  # theta_hat <- solve(t(Z_hat) %*% Z_hat) %*% t(Z_hat) %*% Y
  
  # alternative: direct estimation
  theta_hat <- solve(t(Q) %*% Z) %*% t(Q) %*% Y
  
  return(list(
    lambda = theta_hat[1],
    beta = theta_hat[-1]
  ))
}



####################################################
# Simulation Studies
####################################################


# Simulation parameters
n0_vec <- c(7, 10, 20)  # Different grid sizes
n_rep <- 1000           # Number of replications
lambda_true <- 0.4
beta_true <- c(1, 1, 1) # Intercept, beta1, beta2

# Storage for summary results
simulation_table <- data.frame()

set.seed(2025) # Reproducibility

for(n0 in n0_vec){
  n <- n0^2
  cat(sprintf("Processing sample size: n = %d (grid %dx%d)\n", n, n0, n0))
  
  # 1. Spatial Weight Matrix Construction (Rook Contiguity)
  # -----------------------------------------------------
  coords <- expand.grid(row=1:n0, col=1:n0)
  W0 <- matrix(0, n, n)
  
  # Efficiently build adjacency
  for(i in 1:n){
    for(j in 1:n){
      if(i != j){
        dist <- abs(coords[i,1] - coords[j,1]) + abs(coords[i,2] - coords[j,2])
        if(dist == 1){
          W0[i,j] <- 1
        }
      }
    }
  }
  
  # Row normalization
  rs <- rowSums(W0)
  W <- W0
  if(any(rs > 0)){
    W[rs>0, ] <- W0[rs>0, ] / rs[rs>0]
  }
  
  # Pre-calculate (I - lambda*W)^-1 for DGP
  S_inv_true <- solve(diag(n) - lambda_true * W)
  
  # 2. Monte Carlo Loop
  # -------------------
  # Store estimates: Columns = [lambda, beta0, beta1, beta2]
  res_2sls1 <- matrix(NA, n_rep, 4)
  res_2sls2 <- matrix(NA, n_rep, 4)
  res_biv   <- matrix(NA, n_rep, 4)
  
  for(r in 1:n_rep){
    # Generate X: 1st col is 1, rest are N(0,1)
    X <- cbind(1, rnorm(n), rnorm(n))
    
    # Generate error term V: t-distribution with df=5
    V <- rt(n, df = 5)
    
    # Generate dependent variable Y
    # Y = (I - lambda W)^-1 (X beta + V)
    Y <- S_inv_true %*% (X %*% beta_true + V)
    
    # Method 1: 2SLS (X, WX*)
    est1 <- SAR_2SLS1(Y, X, W)
    res_2sls1[r, ] <- c(est1$lambda, est1$beta)
    
    # Method 2: 2SLS (X, WX*, W^2X*)
    est2 <- SAR_2SLS2(Y, X, W)
    res_2sls2[r, ] <- c(est2$lambda, est2$beta)
    
    # Method 3: Best IV
    est3 <- SAR_BIV(Y, X, W)
    res_biv[r, ] <- c(est3$lambda, est3$beta)
  }
  
  # 3. Calculate Bias and RMSE
  # --------------------------
  # True parameter vector
  theta_true <- c(lambda_true, beta_true)
  
  calc_metrics <- function(est_matrix, true_vec){
    bias <- colMeans(est_matrix) - true_vec
    rmse <- sqrt(colMeans((est_matrix - matrix(true_vec, nrow=nrow(est_matrix), ncol=length(true_vec), byrow=TRUE))^2))
    return(c(bias, rmse))
  }
  
  m1 <- calc_metrics(res_2sls1, theta_true)
  m2 <- calc_metrics(res_2sls2, theta_true)
  m3 <- calc_metrics(res_biv, theta_true)
  
  # Organize into rows for the table
  # Structure: Method, n, Bias_lambda, RMSE_lambda, Bias_beta1, RMSE_beta1 (showing selected params)
  
  # For brevity in the table, let's show stats for lambda and the first slope coefficient (beta_1)
  # Indices in vectors: 
  # Est columns: 1(lambda), 2(beta0), 3(beta1), 4(beta2)
  # Metric vector: 1-4 (Bias), 5-8 (RMSE)
  
  row1 <- data.frame(Method="2SLS1", n=n, 
                     Bias_Lam=m1[1], RMSE_Lam=m1[5], 
                     Bias_B1=m1[2], RMSE_B1=m1[6],
                     Bias_B2=m1[3], RMSE_B2=m1[7], 
                     Bias_B3=m1[4], RMSE_B3=m1[8])
  row2 <- data.frame(Method="2SLS2", n=n, 
                     Bias_Lam=m2[1], RMSE_Lam=m2[5], 
                     Bias_B1=m2[2], RMSE_B1=m2[6],
                     Bias_B2=m2[3], RMSE_B2=m2[7], 
                     Bias_B3=m2[4], RMSE_B3=m2[8])
  row3 <- data.frame(Method="BestIV", n=n, 
                     Bias_Lam=m3[1], RMSE_Lam=m3[5], 
                     Bias_B1=m3[2], RMSE_B1=m3[6],
                     Bias_B2=m3[3], RMSE_B2=m3[7], 
                     Bias_B3=m3[4], RMSE_B3=m3[8])
  
  simulation_table <- rbind(simulation_table, row1, row2, row3)
}

# Print the Final Summary Table
print(simulation_table, digits=4)


####################################################
# Comparison & Analysis
####################################################

```