# iFunMed: Integrative Functional Mediation Analysis of GWAS and eQTL

## Overview

*iFunMed* is a mediation model that utilizes functional annotation data as prior information and builds on summary statistics from GWAS and eQTL  studies. *iFunMed* model capitalizes on the functional annotation information when modeling the probability that a given SNP has a non-zero direct or indirect effect and, as a result, enables identification of SNPs that are associated with phenotypical changes through direct  phenotype-genotype  and/or indirect  phenotype-genotype through gene expression effect. Furthermore, we propose a pipeline to enable selection of the functional annotations for direct and indirect effect models based on enrichment measurements. 

![iFunMed Model Depiction](figures/ModelDepiction.png)


## iFunMed Model Fitting

We provide an example data file (`example_data.RData`) for new users to get familiar with *iFunMed*. It consist of small set of 500 SNPs and contains the input for the model:
- GWAS.summstat: GWAS summary statistics
- eQTL.summstat: eQTL summary statistics
- LD.matrix: LD matrix (500x500)
- annotation.matrix: Functional annotation matrix with 5 different binary annotations (500x5)

`iFunMed_base.r` has the necessary elemets to run the model and consists, mainly, of four functions: 
- `vemMedSS`: Fit the direct effect part of the mediation model, adjusted by the mediator (![EqnMedModel](http://latex.codecogs.com/gif.latex?Z_Y%3D%5CSigma%20%5Cbeta%20&plus;%20Z_G%20%5Cgamma&plus;%20%5Cepsilon)).
- `vemDirectSS`: Fit the gene effect part of the mediation model (![EqnGeneModel](http://latex.codecogs.com/gif.latex?Z_G%3D%5CSigma%20B%20&plus;%20%5Ceta)).
- `processFunMed`: Summarize output. 
- `annoEnrich`: Calculates enrichment values for the annotation matrix using the *iFunMed* fit without annotation.


Once you read the R code and load the data into your R session, you can fit *iFunMed*. You can run the model with annotation directly, or following the annotation selection pipeline.

```
source('iFunMed_base.r')
load('example_data.RData')
```

### 1. Annotation Model

Let's say that out of the 5 annotations, we are only interested in the fifth annotation (`A5`). `anno.iFunMed` will be a 500x2 annotation matrix where the first column represents the intercept (all ones) and the second column if the `A5` binary annotation information. 
Then, we fit the model in a two step fashion with the functions `vemDirectSS` and `vemMedSS`.

```
anno.iFunMed <- cbind(1, annotation.matrix[, 'A5'])
vem.GEM.anno <- vemDirectSS(LD.matrix, GEMsummstats = eQTL.summstat, anno = anno.iFunMed)
vem.DEM.anno <- vemMedSS(LD.matrix, DEMsummstats = GWAS.summstat, GEMsummstats = eQTL.summstat, anno = anno.iFunMed)
```

`vem.GEM.anno` and `vem.DEM.anno` are objects that represent the direct and indirect effect information. We can use the `processFunMed` function to summarize the main outputs:
```
iFunMedNull.output <- processFunMed(GEMoutput = vem.GEM.anno, DEMoutput = vem.DEM.anno)
```

### 2. Annotation Selection Pipeline

#### 2.1 Fit Model Without Annotation (Null Model)

Similarly as in Section 1., we fit the model in a two step fashion with the functions `vemDirectSS` and `vemMedSS`. Since we are fitting the null model (without annotation), we set `anno = NULL` in both of them, which is also the default.
```
vem.GEM.null <- vemDirectSS(LD.matrix, GEMsummstats = eQTL.summstat, anno = NULL)
vem.DEM.null <- vemMedSS(LD.matrix, DEMsummstats = GWAS.summstat, GEMsummstats = eQTL.summstat, anno = NULL)
```

We use `processFunMed` function to summarize the main outputs:
```
iFunMedNull.output <- processFunMed(GEMoutput = vem.GEM.null, DEMoutput = vem.DEM.null)
```

#### 2.2 Measure Annotation Enrichment 

Annotation enrichment values are calculated with the `annoEnrich` function. Based on the the direct and indirect effect posterior probabilities from the null model object (`iFunMedNull.output`) it calculates the average posterior probability of inclusion of the SNPs with the annotation (`avePP`: ![EqnAnnoB](http://latex.codecogs.com/gif.latex?ave%28%5Chat%7Bs%7D_%7BB%2C%20k%7D%29) and ![EqnAnnoBeta](http://latex.codecogs.com/gif.latex?ave%28%5Chat%7Bs%7D_%7B%5Cbeta%2C%20k%7D%29) in the manuscript) and obtains the measurement of enrichment for each annotation (`enrichment`: ![EqnAnnoEnrichB](http://latex.codecogs.com/gif.latex?%5Chat%7Bp%7D_%7BB%7D) and ![EqnAnnoEnrichBeta](http://latex.codecogs.com/gif.latex?%5Chat%7Bp%7D_%7B%5Cbeta%7D) in the manuscript).

The `parallel` library is required for utilizing `annoEnrich`.
```
library(parallel)
annoSelection <- annoEnrich(FunMedNull = iFunMedNull.output, annoMtx = annotation.matrix, Nperm = 10000, cores = 20)
> annoSelection
$avePP
$avePP$GEM
          A1           A2           A3           A4           A5 
1.023347e-06 9.965326e-03 1.350106e-02 1.023347e-06 9.227229e-03 

$avePP$DEM
          A1           A2           A3           A4           A5 
1.176456e-02 4.740930e-11 6.992922e-03 1.923054e-02 9.259147e-03 


$enrichment
$enrichment$GEM
    A1     A2     A3     A4     A5 
0.9971 0.2232 0.0770 0.9826 0.2499 

$enrichment$DEM
    A1     A2     A3     A4     A5 
0.3114 0.9972 0.4824 0.1963 0.3829 
```

Since the enrichment is calculated based on permutation, these values may vary slightly among different runs. 

#### 2.3 Fit Model With Most Enriched Annotations 

Based on the previous values, we can fit *iFunMed* with the most enriched annotation for each model and add the intercept to `A3` and `A4`. Alternatively, multiple annotations can be used.
```
annomtx.GEM <- cbind(1, annotation.matrix[, "A3"])
annomtx.DEM <- cbind(1, annotation.matrix[, "A4"])

vem.GEM.anno <- vemDirectSS(LD.matrix, GEMsummstats = eQTL.summstat, anno = annomtx.GEM )
vem.DEM.anno <- vemMedSS(LD.matrix, DEMsummstats = GWAS.summstat, GEMsummstats = eQTL.summstat, anno = annomtx.DEM)
```

Similarly, we use `processFunMed` to summarize results:
```
processFunMed(GEMoutput = vem.GEM.anno, DEMoutput = vem.DEM.anno)
```

### 3. `processFunMed` Output and Model Parameters

From the fitting without and with annotation, `processFunMed` will summarize the information from `vemDirectSS` and `vemMedSS`.
`iFunMedNull.output` and `iFunMedAnno.output` are lists with the following information:

- Convergency: Number of iterations and convergency status for direct and indirect efect models.
- Parameters: Direct and indirect effect estimated model parameters.
- PostProb: Posterior Probability of  inclusion (non-zero effect size) and FDR-corrected values for direct and indirect efect models.

```
> str(iFunMedAnno.output)
List of 3
 $ Convergency:List of 2
  ..$ GEM:List of 2
  .. ..$ niter    : num 19
  .. ..$ converged: logi TRUE
  ..$ DEM:List of 2
  .. ..$ niter    : num 14
  .. ..$ converged: logi TRUE
 $ Parameters :List of 2
  ..$ GEM:List of 3
  .. ..$ gammaB: num [1:2] -5.466 0.791
  .. ..$ varEta: num 1.05
  .. ..$ nuB   : num 17.3
  ..$ DEM:List of 4
  .. ..$ gammaBeta : num [1:2] -5.97 1.93
  .. ..$ varEpsilon: num 0.886
  .. ..$ nuBeta    : num 26.9
  .. ..$ gamma     : num 0.0139
 $ PostProb   :List of 2
  ..$ GEM:List of 2
  .. ..$ PP    : Named num [1:500] 1.67e-06 7.56e-07 7.56e-07 1.67e-06 7.56e-07 ...
  .. .. ..- attr(*, "names")= chr [1:500] "SNP1" "SNP2" "SNP3" "SNP4" ...
  .. ..$ FDR.PP: Named num [1:500] 0.823 0.994 0.994 0.965 0.988 ...
  .. .. ..- attr(*, "names")= chr [1:500] "SNP1" "SNP2" "SNP3" "SNP4" ...
  ..$ DEM:List of 2
  .. ..$ PP    : Named num [1:500] 2.52e-08 3.68e-09 3.68e-09 2.52e-08 3.68e-09 ...
  .. .. ..- attr(*, "names")= chr [1:500] "SNP1" "SNP2" "SNP3" "SNP4" ...
  .. ..$ FDR.PP: Named num [1:500] 0.92 0.987 0.993 0.894 0.985 ...
  .. .. ..- attr(*, "names")= chr [1:500] "SNP1" "SNP2" "SNP3" "SNP4" ...
```






