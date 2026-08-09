Εργασία 3
================
Παναγίωτης Δημητρίου
2026-04-28

## ΑΣΚΗΣΗ 1

Γνωρίζουμε από ΙΝΜΑ ότι ο εκτιμητης Monte-Carlo
$Z_n^{MC} \to \mathbb{E}[h(X)]$ όπου Xj ακολουθούν την $F$. Σε αυτή την
περίπτωση έχουμε ότι $Xj\sim Exp(2)$ και $h(x) = \sin(x) \log^3(x)$.
Αρχικά εκτιμήσαμε την διασπορά με την εντολή της R var() για τυχαίο
δείγμα μεγέθους m=100000 από την Exp(2). Έπειτα, επειδή
$MSE(z_n^{MC}) = \frac{Var(h(X))}{n}$, και θέλουμε επιθυμητό σφάλμα της
τάξης του $\epsilon > 0$, τότε
$$\frac{Var(h(X))}{n} < \epsilon^2 \implies n \ge \frac{Var(h(X))}{\epsilon^2} \quad (1)$$
Ισχύει ότι
$$RMSE(z_n^{MC}) = \sqrt{MSE(z_n^{MC})} = \frac{\sqrt{Var(h(X))}}{\sqrt{n}} \sim \text{σφάλμα}$$
και έτσι θα βρούμε την κατάλληλη τιμή n που είναι το μέγεθος του
δείγματος για να προσομειώσουμε απο την Exp(2) έτσι ώστε να βγάζουμε ένα
Monte-Carlo.

``` r
mctheta=function(N,eps){
  
  if(N<=0 || N%%1!=0){
    stop("Το N θα πρέπει να είναι θετικός ακέραιος")
  }
  if(eps<=0){
    stop("Το eps θα πρέπει να είναι θετικός")
  }
  x1=rexp(100000,2)
  x2=sin(x1)*(log(x1))^3 #εκτίμηση του var(h(X))
  esvar=var(x2)
  
  n=ceiling(esvar/eps) #αφού n>=var(h(X))/ε^2 τότε παίρνω ceiling για τον επόμενο ακέραιο n που ισχυεί για το επιθημητό σφάλμα ε και βρίσκουμε το πλήθος n τετοιό ώστε το RMSE<ε
  sim.MC=c()
  for(i in 1:N){ #εκτελούμε το MC πείραμα Ν φορές
  z=rexp(n,2)
  mc=mean(sin(z)*log(z)^3) #υπολογισμός Monte-Carlo εκτιμητή με Ζ να ακολουθούν την Exp(2) μεγέθους n στο πλήθος, όπως αυτό υπολογίστηκε στο προηγούμενο στάδιο για να έχουμε RMSE<ε
  sim.MC[i]=mc #διάνυσμα από της Ν εκτιμήσεις Monte-Carlo
  }
  
  return(list(vector=sim.MC,sample.size=n))
}
```

Θα εκτελέσουμε το function για σφάλματα ϵ = 0.01,0.001,0.0001 και N =
100 (αριθμός επαναλήψεων πειράματος Monte-Carlo) το οποίο θα μας
επιστρέψει ένα διάνυσμα μεγέθους N=100 με τα δείγματα της εκτιμήτριας
MC. Σε κάθε περίπτωση αναλόγως του ε, το n λαμβάνει διαφορετική τιμή
όπως αναφέραμε.

``` r
set.seed(1)
z1=mctheta(100,0.01)
z2=mctheta(100,0.001)
z3=mctheta(100,0.0001)
```

Γενικά από την $(1)$ αναμένουμε ότι καθώς το ε γίνεται μικρότερο τότε
θέλουμε μεγαλύτερο δείγμα προσομείωσης.

Η Monte-Carlo εκτιμήτρια $z_n^{MC}$ του $\mathbb{E}[h(X)]$ είναι
αμερόληπτη
($\mathbb{E}(z_n^{MC}) = \mathbb{E}[h(X)] \iff bias(z_n^{MC}) = 0$) και
έχει μέσο τετραγωνικό σφάλμα:

$$MSE(z_n^{MC}) = Var(z_n^{MC}) = \frac{1}{n} Var(h(X))$$

**Απόδειξη:**

$$
\begin{aligned}
\mathbb{E}(z_n^{MC}) &= \mathbb{E}\left(\frac{1}{n} \sum_{j=1}^{n} h(X_j)\right) \\
&= \frac{1}{n} \sum_{j=1}^{n} \mathbb{E}(h(X_j)) \stackrel{iid}{=} \frac{1}{n} n \cdot \mathbb{E}[h(X)] = \mathbb{E}[h(X)]
\end{aligned}
$$

Για το MSE (εκτιμητής για σφάλμα στο τετράγωνο / τετράγωνο του
σφάλματος) έχουμε το ακόλουθο:

$$
\begin{aligned}
MSE(z_n^{MC}) &= Var(z_n^{MC}) = Var\left(\frac{1}{n} \sum_{j=1}^{n} h(X_j)\right) \\
&\stackrel{iid}{=} \frac{1}{n^2} \sum_{j=1}^{n} Var(h(X)) = \frac{n}{n^2} Var(h(X)) = \frac{1}{n} Var(h(X))
\end{aligned}
$$ Και απο Κεντρικό Οριακό Θεώρημα ισχύει.
$$Z_n^{MC} \sim N\left(\mathbb{E}[h(X)], \frac{Var(h(X))}{n}\right)$$

``` r
mcn1=z1$vector 
n1=z1$sample.size
n1 #χρειαστίκαμε δείγμα μήκους n=26 για το ε=0.01
```

    ## [1] 26

``` r
mcn2=z2$vector
n2=z2$sample.size
n2
```

    ## [1] 255

``` r
mcn3=z3$vector
n3=z3$sample.size
n3
```

    ## [1] 2547

Πιο κάτω θα κάνουμε τα αντίστοιχα ιστογράμματα για το κάθε διάνυσμα που
βγάλαμε έχοντας τον ίδιο αριθμό Μonte-Carlo πειραμάτων(N=100).

``` r
hist(mcn1, breaks = 15,xlim=c(-0.75,-0.2),main="samples of MC estimator with n1")
```

![](Statistical-Data-Analysis-Project%5BR%5D-3-_files/figure-gfm/unnamed-chunk-4-1.png)<!-- -->

``` r
hist(mcn2, breaks = 15,xlim=c(-0.65,-0.4),main="samples of MC estimator with n2")
```

![](Statistical-Data-Analysis-Project%5BR%5D-3-_files/figure-gfm/unnamed-chunk-4-2.png)<!-- -->

``` r
hist(mcn3, breaks = 15,xlim=c(-0.65,-0.4),main="samples of MC estimator with n3")
```

![](Statistical-Data-Analysis-Project%5BR%5D-3-_files/figure-gfm/unnamed-chunk-4-3.png)<!-- -->

Βλέπουμε ότι όπως περιμέναμε από τη θεωρία, καθώς το n μεγαλώνει τα
δείγματα συγκεντρώνονται γύρω από μια μέση τιμή και η διακύμανση
ελαττώνεται κατά $1/n$. Με την αλλαγή του κάθε δείγματος που αντιστοιχεί
σε διαφορετικό n το οποίο μεγαλώνει έχουμε “καλύτερο” ιστόγραμμα αφού η
κατανομή της εκτιμήτριας συγκλίνει προς την κανονική όπως αναμέναμε.

``` r
par(mfrow=c(1,2))
boxplot(mcn1,main="Boxplot mcn1")
boxplot(mcn2,main="Boxplot mcn2")
```

![](Statistical-Data-Analysis-Project%5BR%5D-3-_files/figure-gfm/unnamed-chunk-5-1.png)<!-- -->

``` r
par(mfrow=c(1,1))
boxplot(mcn3,ylim=c(-0.55,-0.48),main="Boxplot mcn3")
```

![](Statistical-Data-Analysis-Project%5BR%5D-3-_files/figure-gfm/unnamed-chunk-5-2.png)<!-- -->

Γενικά καθώς μεγάλωνε το n τα boxplot σε κάθε διαφορετικό δείγμα
mcn1(n=26), mcn2(n=254), mcn3(n=2545) φαίνοντε να μοιάζουν ακόμη
περισσότερο στο θηκόγραμμα της κανονικής.Για παράδειγμα το lower whisker
για τον mcn1 είναι πολύ μεγαλύτερο ενώ στο mcn3 έχει ίδιο μήκος με το
upper whisker. Επίσης ισχύει και για το mcn2 (έχει και outlier). Δηλαδή
παρουσιάζουν κάποιες αποκλίσεις σε θέμα σχήματος.

``` r
smcn1=scale(mcn1)
ks.test(smcn1,"pnorm")
```

    ## 
    ##  Asymptotic one-sample Kolmogorov-Smirnov test
    ## 
    ## data:  smcn1
    ## D = 0.055383, p-value = 0.9189
    ## alternative hypothesis: two-sided

``` r
smcn2=scale(mcn2)
ks.test(smcn2, 'pnorm')
```

    ## 
    ##  Asymptotic one-sample Kolmogorov-Smirnov test
    ## 
    ## data:  smcn2
    ## D = 0.054664, p-value = 0.9262
    ## alternative hypothesis: two-sided

``` r
smcn3=scale(mcn3)
ks.test(smcn3, 'pnorm')
```

    ## 
    ##  Asymptotic one-sample Kolmogorov-Smirnov test
    ## 
    ## data:  smcn3
    ## D = 0.037347, p-value = 0.999
    ## alternative hypothesis: two-sided

Κανονικοποιήσαμε το διάνυσμα μας και με την Kolmogorov-Smirnov
παρατηρούμε ότι το smcn1(n=26) και το smcn2(n=255) ενώ είχαν εξαιρετικό
p-value 0.9189 και 0.9262 αντίστοιχα, που σημαίνει ότι ακολουθούν την
κανονική το smcn3 είναι ακόμη πιο τέλειο.

``` r
qqnorm(smcn1)
qqline(smcn1)
```

![](Statistical-Data-Analysis-Project%5BR%5D-3-_files/figure-gfm/unnamed-chunk-7-1.png)<!-- -->

``` r
qqnorm(smcn2)
qqline(smcn2)
```

![](Statistical-Data-Analysis-Project%5BR%5D-3-_files/figure-gfm/unnamed-chunk-7-2.png)<!-- -->

``` r
qqnorm(smcn3)
qqline(smcn3)
```

![](Statistical-Data-Analysis-Project%5BR%5D-3-_files/figure-gfm/unnamed-chunk-7-3.png)<!-- -->

Επίσης απο το qqplots βλέπουμε ότι στις πρώτες δύο περιπτώσεις είχαμε τα
σημεία μας είχαν πιο μεγάλη απόσταση απο την τελευταία περίπτωση.

``` r
var(mcn1)/var(mcn2)
```

    ## [1] 6.702366

``` r
var(mcn1)/var(mcn3)
```

    ## [1] 68.75365

``` r
var(mcn2)/var(mcn3)
```

    ## [1] 10.25812

Από όλα αυτά επιβεβαιώνουμε την σύγκλιση της Monte-Carlo εκτιμήτριας
προς την κανονική κάθως το n μεγαλώνει αλλά και την μείωση της διασποράς
απο το προηγούμενο πηλίκο που το αποδεικνύει.

## ΑΣΚΗΣΗ 2

Θέλετε να υπολογίσετε τη μέση τιμή $E[h(X)]$ για την τυχαία μεταβλητή
X∼N(0,1) με συνάρτηση
$$h(x) = e^{-\frac{(x-3)^2}{2}} + e^{-\frac{(x-6)^2}{2}}$$

(β) Γράψτε ένα πρόγραμμα στην R για να την εκτιμήσετε χρησιμοποιώντας
δείγμα μεγέθους $n = 10^3$ και υπολογίστε το σφάλμα.

``` r
set.seed(123)
n = 1000
x = rnorm(n,0,1) #δημιουργούμε ένα δείγμα από την κανονική κατανομή

h =function(x){
  exp(-(x-3)^2 / 2) + exp(-(x-6)^2 / 2)
}

h_values = h(x) #υπολογισμός των h(x)
mc_estimate = mean(h_values) #υπολογισμός Monte Carlo εκτίμησης της θεωρητικής ποσότητας μ=Ε[h(x)].

#Υπολογισμός του τυπικού σφάλματος (Standard Error)
#Το σφάλμα της Monte Carlo εκτίμησης είναι s / sqrt(n)
standard_error = sd(h_values) / sqrt(n)

# Εκτύπωση αποτελεσμάτων
cat("Εκτίμηση Monte Carlo:", mc_estimate, "\n")
```

    ## Εκτίμηση Monte Carlo: 0.07699182

``` r
cat("Τυπικό Σφάλμα:", standard_error, "\n")
```

    ## Τυπικό Σφάλμα: 0.005086384

Από το ερώτημα (α), βρήκαμε ότι
$$ E[h(X)] =\frac{1}{\sqrt{2}}\left(e^{-9/4}+e^{-9}\right)\approx 0.07462 $$,
επομένως έχουμε ότι η εκτίμηση μέσω της Monte Carlo εκτιμήτριας είναι
αρκετά καλή, αφού $$ \hat{\mu}=0.07699 $$ με τυπικό σφάλμα (standard
error) $$ SE=0.00509 $$ και
$$ \left|\hat{\mu}-E[h(X)]\right| = |0.07699-0.07462| = 0.00237 $$ που
είναι μικρότερο από το standard error, άρα το αποτέλεσμα θεωρείται
ικανοποιητικό.

(γ) Υπολογείστε τη μέση τιμή, γράφοντας ένα πρόγραμμα στην R με τη χρήση
του Importance Sampling και χρησιμοποιώντας ως proposal g την N(2.5,2)
μία μίξη κανονικών $1/2*N(1.5,0.5) + 1/2*N(3,0.5)$ καθώς και την
ομοιόμορφη U(−8,−1) και δείγμα μεγέθους $N = 10^3$.

Ποια proposal έχει την καλύτερη απόδοση;

Θα ορίσουμε πρώτα τις συναρτήσεις h(x) και f(x).

``` r
set.seed(42)
N=1000
# ορίζω τα h(x) και f(x)
h = function(x){
  exp(-(x - 3)^2 / 2) + exp(-(x - 6)^2 / 2)}
f =function(x){
  dnorm(x,0,1)}
```

Γενικά, ισχύει
$$(1) RMSE(z_n^{IS}) = \sqrt{Var(z_n^{IS})}=\frac{\sqrt{Var\left(h(X_j)  \frac{f(X_j)}{g(X_j)}\right)}}{\sqrt{n}}$$

Το οποίο είναι το σφάλμα

Άρα η Proposal που δίνει την μικρότερη διακύμανση τότε θα έχουμε το
μικρότερο σφάλμα(το n θα είναι 1000) άρα και καλύτερη απόδοση.

# Proposal 1: N(2.5, 2)

``` r
imps1=function(n){

prop1_samples = rnorm(n, 2.5, sqrt(2))
weights1 = (h(prop1_samples) * f(prop1_samples)) / dnorm(prop1_samples, 2.5,sqrt(2))     

return(mean(weights1))
}
```

``` r
v1=c()
for(i in 1:1000){
  v1[i]=imps1(1000)
  
}
```

``` r
estim.is1=mean(v1)
estim.is1 #εκτίμιση μέσης τιμής με την proposal 1
```

    ## [1] 0.07449862

``` r
var1=var(v1)
var1
```

    ## [1] 5.41721e-06

# Proposal 2: Mixture 0.5*N(1.5, 0.5) + 0.5*N(3, 0.5)

Επιλογή από ποια κανονική θα πάρουμε δείγμα:

``` r
imps2=function(n){
prop2_samples = 1/2*rnorm(n, 1.5, sqrt(0.5))+1/2*rnorm(n, 3, sqrt(0.5))
g2_pdf = 0.5 * dnorm(prop2_samples, 1.5, sqrt(0.5)) + 0.5 * dnorm(prop2_samples, 3, sqrt(0.5))
weights2 = (h(prop2_samples) * f(prop2_samples)) / g2_pdf
return(mean(weights2))
}
```

``` r
v2=c()
for(i in 1:1000){
  v2[i]=imps2(1000)
  
}
```

``` r
estim.is2=mean(v2)
estim.is2
```

    ## [1] 0.07459923

``` r
var2=var(v2)
var2
```

    ## [1] 1.563113e-06

# Proposal 3: U(-8, -1)

Η ομοιόμορφη $U(-8,-1)$ δεν αποτελεί κατάλληλη proposal κατανομή, διότι
δεν καλύπτει την περιοχή όπου η $h(x)f(x)$ παίρνει σημαντικές τιμές, άρα
παραβιάζει τη βασική προϋπόθεση του Importance Sampling.

Αυτό συμβαίνει επειδή f(x)\>0 για κάθε x∈R και h(x)\>0 για κάθε x∈R,
επομένως h(x)f(x)\>0 για κάθε x∈R, με αποτέλεσμα να πρέπει η proposal
κατανομή μας να παίρνει θετικές τιμές σχεδόν παντού στα σημεία αυτά,
πράγμα που δεν ισχύει σε αυτή την περίπτωση.

Συγκεκριμένα, μηδενίζεται σε περιοχές (όπως γύρω από το 3 και το 6) όπου
το γινόμενο $h(x)f(x)$ παίρνει σημαντικές θετικές τιμές, παραβιάζοντας
έτσι τη βασική θεωρητική προϋπόθεση του Importance Sampling ($g(x) > 0$
όταν $h(x)f(x) \neq 0$).

Ας το επιβεβαιώσουμε και μέσω των πιο κάτω αποτελεσμάτων:

``` r
imps3=function(n){
prop3_samples = runif(N,-8,-1)
g3_pdf = dunif(prop3_samples,-8,-1)
weights2 = (h(prop3_samples) * f(prop3_samples)) / g3_pdf
return(mean(weights2))
}
```

``` r
v3=c()
for(i in 1:1000){
  v3[i]=imps3(1000)
  
}
```

``` r
estim.is3=mean(v3)
estim.is3
```

    ## [1] 1.518462e-05

``` r
var3=var(v3)
var3
```

    ## [1] 4.276662e-12

Η ακαταλληλότητα της τρίτης proposal κατανομής $g$ επιβεβαιώνεται από τη
μέση τιμή που είναι πολύ μακριά από την θεωρητική μέση τιμή και την
διακύμανση που είναι σχεδόν 0, γεγονός που οδηγεί σε λάθος εκτιμητή.
Επειδή η $h(x)$ στο διάστημα $[-8, -1]$ είναι πρακτικά 0, όλα τα
δείγματα έδωσαν αποτέλεσμα 0.

Συνεπώς, η $U(-8,-1)$ παρουσιάζει κακή απόδοση και αποτέλει λανθασμένη
επιλογή.

Πιο κάτω μπορούμε να δούμε πιο καθαρά τα αποτελέσματα που προκύπτουν από
τις τρεις proposal g μέσω της χρήσης του Importance Sampling:

``` r
results=data.frame(
  Proposal = c("N(2.5, 2)", "Mixture Normal","U(-8,-1)"),
  Variance_of_Est=c(var1,var2,var3)
) 

print(results)
```

    ##         Proposal Variance_of_Est
    ## 1      N(2.5, 2)    5.417210e-06
    ## 2 Mixture Normal    1.563113e-06
    ## 3       U(-8,-1)    4.276662e-12

``` r
var1/var2
```

    ## [1] 3.465655

Η διακύμανση var2 που αντιστοιχεί στην Mixture Normal είναι 3.5 φορές
μικρότερη από την var1 της N(2.5,2).

Eπομένως από την σχέση (1) έχουμε ότι: Από τον πίνακα των διακυμάνσεων
φαίνεται ότι η καλύτερη απόδοση επιτυγχάνεται με την Mixture Normal
κατανομή, καθώς παρουσιάζει το μικρότερο σφάλμα (variance) ανάμεσα στις
έγκυρες επιλογές (1.56e-06). Η ομοιόμορφη κατανομή $U(-8,-1)$ δίνει
μηδενική διακύμανση (4.27e-12), διότι παίρνει δειγματα σε περιοχή όπου η
συνάρτηση μηδέν, οδηγώντας σε έναν λανθασμένο εκτιμητή.
