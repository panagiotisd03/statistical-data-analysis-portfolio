# Applied Statistical Modeling & Analysis in R

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
*  **Repository Files**
| File | Description |
| :--- | :--- |
| [View pdf](<./project-1/Statistical Data Analysis Project[R](1).pdf>) | The primary pdf with R code|
---

###  Time Series Simulation, Inverse Transform Sampling ,Acceptance-Rejection Sampling, Dataset Analysis
* **AR(1) Time Series Process:**
  * Implementation of a custom simulator for $X_t = \mu + \rho X_{t-1} + \epsilon_t$.
  * Simulate the time series across varying parameter configurations, testing the effects of different initial states , mean shifts, autocorrelation coefficients, and
  * Assessment of stationarity, autocorrelation at lag-1, and normality verification (Kolmogorov-Smirnov normality test, Q-Q plots,kernel,ecdf).
* **Data Analysis of `iris` dataset:**
  * Summary metrics (mean, variance, IQR), outlier detection, and goodness-of-fit distribution testing.
  * **Random Variate Generation:**  Generated random samples using the Inverse Transform Sampling method.
  * **Acceptance-Rejection Sampling:**  Implemented the Acceptance-Rejection method using heavy-tailed proposal distributions to simulate exponential random variables and analyzed sampling efficiency.
**Repository Files**
| File | Description |
| :--- | :--- |
[View RMarkdown Code](<./project-2/Statistical Data Analysis Project[R](2).Rmd>)| The primary RMarkdown source code |
| [View](<./project-2/Statistical-Data-Analysis-Project[R]-2-.md>) | The final form|
---

###  Monte-Carlo Simulation, Importance Sampling
* Implemented Monte Carlo integration
* Evaluated Importance Sampling across multiple proposal distributions (Gaussian, Gaussian Mixture, Uniform) for variance reduction and convergence analysis.
*Compared Standard MC against 3 proposal distributions were the Mixture Normal proposal yielded the lowest variance
**Repository Files**
| File | Description |
| :--- | :--- |
| [View RMarkdown Code](<./project-3/Statistical Data Analysis Project[R](3).Rmd>) | The primary RMarkdown source code |
| [View](<./project-3/Statistical-Data-Analysis-Project[R]-3-.md>)  | The final form|
---

##  Languages & Tools

* **Language:** R
* **Libraries:** `moments`, `nortest`, `rmutil`, `ggplot2`, `stats`, `knitr`
* **Workflow:** RMarkdown, RStudio

---
