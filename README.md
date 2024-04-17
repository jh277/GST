**GST package**

**Package description**:

Determining whether one treatment is preferred over others is a major goal of many clinical studies 
but can be complicated by the situation when no single outcome is sufficient to make the judgement. 

This package presents a useful global statistics test (GST) technique and the corresponding global treatment effect (GTE) measure for assessing treatment's global preference when multiple outcomes are evaluated together.

The package contains 6 functions:

*\. GST.parameter() provides unbiased estimators of GST parameters for public sharing.
*\.  GST() performs one or two-sided non-parametric global statistical test.
*\.  theta() calculates global treatment effect.
*\.  N.GSTwPilotData() calculates sample size needed in each group when there is no pilot data.
*\.  N.GSTnoData() calculates sample size when some pilot data are available.
*\.  GST.simu() performs simulation using GST


**Installation and use**:

1. In your R console, run "library(devtools)". If you have not downloaded devtools package, please run "install.packages("devtools")" first.

2. In your R console, run "devtools::install_github("jh277/GST")" to install this R package.

3. run "library(GST)" in your console to make it available in your current R session.

