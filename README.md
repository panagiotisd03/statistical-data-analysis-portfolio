# Data Analysis

This repository contains my coursework, projects, and simulations in a course called Statistical Data Analysis.
*University of Cyprus — Department of Mathematics and Statistics* 
---
 R Fundamentals, Data Analysis, Visualization 
---
* **R Fundamentals:**
Functions,loops, matrices
* **Data Analysis:**
  * Samples, visualization of the sample[hist(), boxplot(),ECDF plots()],mean(),median() and conclusions [skewness, quantiles]
* **Laplace Distribution Analysis:**
  * **Simulation:** Random sample generation from the Laplace distribution.
  * **Goodness-of-Fit Testing:** Evaluated goodness-of-fit via the Kolmogorov-Smirnov (K-S) test, Gaussian Kernel Density Estimation (KDE), and Q-Q plots against theoretical quantiles.
  * **Theoretical vs Empirical Moments:** Analytical integration to derive theoretical expectation followed by comparison with sample statistics.[hist(),density(),ecdf()]
---

###  Time Series Simulation, EDA & Acceptance-Rejection Sampling
* **AR(1) Time Series Process:**
  * Implementation of a custom simulator for $X_t = \mu + \rho X_{t-1} + \epsilon_t$.
  * Assessment of stationarity, autocorrelation at lag-1, and normality verification (Shapiro-Wilk test, Q-Q plots).
* **Statistical Profiling (`iris` dataset):**
  * Summary metrics (mean, variance, IQR), outlier detection, and goodness-of-fit distribution testing.
* **Non-Uniform Random Variate Generation:**
  * **Inverse Transform Method:** Analytical inversion of $F(x) = \frac{x}{x+1}$ ($x > 0$), handling infinite moments ($\infty$), and ECDF verification.
  * **Acceptance-Rejection Sampling:** Mathematical proof and implementation to generate $\text{Exp}(\lambda)$ variates using an auxiliary envelope density function. Evaluated rejection rates and performed goodness-of-fit tests on trial counts against the Geometric distribution.

---

### 🔹 Assignment 3: Advanced Topics

---

## 🛠️ Languages & Tools

* **Language:** R
* **Libraries:** `moments`, `nortest`, `rmutil`, `ggplot2`, `stats`, `knitr`
* **Workflow:** RMarkdown, VS Code / Codespaces, Git

---
