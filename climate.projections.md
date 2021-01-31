Goal
----

To outline three methods to consider how one could simulate future
states of an environmental variable (E) based on past observations or
models. The goal of such work would be to simulate future E to drive say
a population dynamics model that has a functional E dependence to
account for a mean change in E owing to climate change.

Outline
-------

Three methods:

1.  Non-parametric resampling

-   blocked sampling
-   quantile sampling

1.  Parametric resampling

-   mean and variance shifted

1.  Down-scaled climate projections

-   an example using the Climate Atlas of Canada

First, a data set
-----------------

We will use the R package gslea to get data for the Gulf of St. Lawrence

    devtools::install_github("duplisea/gslea")
    library(gslea)

This package gathers together lots of different data for broad ecosystem
approach regions of the Gulf of St. Lawrence.

<img src="https://github.com/duplisea/gslea/blob/master/README_files/figure-markdown_strict/gslmap.plain-1.png" style="width:100.0%" />

This package has search functions that we can use to find a variable.
Here we search for “temp” as part of temperature and then we select the
variable T250 (temperature at 250m) and we create a data.table with that
variable for all years available for ecosystem approach region 2
(central Gulf). We also use the search function to create a text string
describing the variable that can be used for axis label and titles:

    find.vars.f("temp")

    ##   [1] "Ann.coolingdd.high.RCP45"    "Ann.coolingdd.high.RCP85"   
    ##   [3] "Ann.coolingdd.low.RCP45"     "Ann.coolingdd.low.RCP85"    
    ##   [5] "Ann.coolingdd.med.RCP45"     "Ann.coolingdd.med.RCP85"    
    ##   [7] "Ann.distr.coolingdd.mean.45" "Ann.distr.coolingdd.mean.85"
    ##   [9] "Ann.distr.coolingdd.sd.45"   "Ann.distr.coolingdd.sd.85"  
    ##  [11] "Ann.distr.mean.fdd.mean.45"  "Ann.distr.mean.fdd.mean.85" 
    ##  [13] "Ann.distr.mean.fdd.sd.45"    "Ann.distr.mean.fdd.sd.85"   
    ##  [15] "Ann.distr.minus15.mean.45"   "Ann.distr.minus15.mean.85"  
    ##  [17] "Ann.distr.minus15.sd.45"     "Ann.distr.minus15.sd.85"    
    ##  [19] "Ann.mean.T.high.RCP45"       "Ann.mean.T.high.RCP85"      
    ##  [21] "Ann.mean.T.low.RCP45"        "Ann.mean.T.low.RCP85"       
    ##  [23] "Ann.mean.T.med.RCP45"        "Ann.mean.T.med.RCP85"       
    ##  [25] "Ann.mean.fdd.high.RCP45"     "Ann.mean.fdd.high.RCP85"    
    ##  [27] "Ann.mean.fdd.low.RCP45"      "Ann.mean.fdd.low.RCP85"     
    ##  [29] "Ann.mean.fdd.med.RCP45"      "Ann.mean.fdd.med.RCP85"     
    ##  [31] "Ann.minus15.high.RCP45"      "Ann.minus15.high.RCP85"     
    ##  [33] "Ann.minus15.low.RCP45"       "Ann.minus15.low.RCP85"      
    ##  [35] "Ann.minus15.med.RCP45"       "Ann.minus15.med.RCP85"      
    ##  [37] "AnnTdistr.mean.T.mean.45"    "AnnTdistr.mean.T.mean.85"   
    ##  [39] "AnnTdistr.mean.T.sd.45"      "AnnTdistr.mean.T.sd.85"     
    ##  [41] "Fal.mean.T.high.RCP45"       "Fal.mean.T.high.RCP85"      
    ##  [43] "Fal.mean.T.low.RCP45"        "Fal.mean.T.low.RCP85"       
    ##  [45] "Fal.mean.T.med.RCP45"        "Fal.mean.T.med.RCP85"       
    ##  [47] "FalTdistr.mean.T.mean.45"    "FalTdistr.mean.T.mean.85"   
    ##  [49] "FalTdistr.mean.T.sd.45"      "FalTdistr.mean.T.sd.85"     
    ##  [51] "J.GSNW.Q1"                   "J.GSNW.Q2"                  
    ##  [53] "J.GSNW.Q3"                   "J.GSNW.Q4"                  
    ##  [55] "SST"                         "SST.anomaly"                
    ##  [57] "SST.month10"                 "SST.month11"                
    ##  [59] "SST.month5"                  "SST.month6"                 
    ##  [61] "SST.month7"                  "SST.month8"                 
    ##  [63] "SST.month9"                  "Spr.mean.T.high.RCP45"      
    ##  [65] "Spr.mean.T.high.RCP85"       "Spr.mean.T.low.RCP45"       
    ##  [67] "Spr.mean.T.low.RCP85"        "Spr.mean.T.med.RCP45"       
    ##  [69] "Spr.mean.T.med.RCP85"        "SprTdistr.mean.T.mean.45"   
    ##  [71] "SprTdistr.mean.T.mean.85"    "SprTdistr.mean.T.sd.45"     
    ##  [73] "SprTdistr.mean.T.sd.85"      "Sum.mean.T.high.RCP45"      
    ##  [75] "Sum.mean.T.high.RCP85"       "Sum.mean.T.low.RCP45"       
    ##  [77] "Sum.mean.T.low.RCP85"        "Sum.mean.T.med.RCP45"       
    ##  [79] "Sum.mean.T.med.RCP85"        "SumTdistr.mean.T.mean.45"   
    ##  [81] "SumTdistr.mean.T.mean.85"    "SumTdistr.mean.T.sd.45"     
    ##  [83] "SumTdistr.mean.T.sd.85"      "T.deep"                     
    ##  [85] "T.shallow"                   "T150"                       
    ##  [87] "T200"                        "T250"                       
    ##  [89] "T300"                        "Tmax200.400"                
    ##  [91] "Win.mean.T.high.RCP45"       "Win.mean.T.high.RCP85"      
    ##  [93] "Win.mean.T.low.RCP45"        "Win.mean.T.low.RCP85"       
    ##  [95] "Win.mean.T.med.RCP45"        "Win.mean.T.med.RCP85"       
    ##  [97] "WinTdistr.mean.T.mean.45"    "WinTdistr.mean.T.mean.85"   
    ##  [99] "WinTdistr.mean.T.sd.45"      "WinTdistr.mean.T.sd.85"     
    ## [101] "H.NAO"

    var="T250"
    vardesc=find.vars.f(var,T)$description
    varval= EA.query.f(var,EARs=2, years=1800:2020)

Now let’s plot the data and look at how it is distributed

    par(mfcol=c(2,1),mar=c(5,5,.3,.3))
    plot(varval$year, varval$value, type= "b", pch=20, xlab="Year", ylab=vardesc)
    plot(density(varval$value),ylab="Density",xlab=vardesc, main="")
    legend("topright",bty="n",cex=0.7, legend=paste("Shapiro-Wilks =",round(shapiro.test(varval$value)$p.value,3)))

![](climate.projections_files/figure-markdown_strict/gsl.data.dist-1.png)

The data looks pretty well normally distributed and you can see that
since about 2010, the temperature is a steady rise to be at peak levels
in the most recent years.

Non-parametric sampling (empirical)
-----------------------------------

So using this dataset we can develop a scheme for sampling from the
observations either in blocks or from quantiles of the distribution
which could be hypothesised to represent a distribution for future
climate.

So as a first example, say one wants to create future climate projection
based on only the last 10 years of data and say 1000 samples for a
projection

    plot(varval$year, varval$value, type= "b", pch=20, xlab="Year", ylab=vardesc)
    Tsamps= sample(tail(varval$value,10),size=1000,replace=T)
    lines(tail(varval$year,10),rep(median(Tsamps),10),col="red")
    lines(tail(varval$year,10),rep(quantile(Tsamps,.75),10),col="red",lty=2)
    lines(tail(varval$year,10),rep(quantile(Tsamps,.25),10),col="red",lty=2)

![](climate.projections_files/figure-markdown_strict/empirical.sampling-1.png)

This shows the median and the 50% confidence interval on the median
based on the resampling of the last 10 years of data. This sample could
form the basis of a plausible near future climate scenario.

The problem with this kind of sampling is that it could be relatively
conservative because it is constrained to observed values. This means
that it may be too conservative to account for plausible future climate
scenarios. It also is constrained to the observed values only while we
might get a smoother and more plausble future probability density
function using interpolation.

parametric resampling
---------------------

The non-parametric approach described above can be ‘smoothed’ by
considering a parametric approach where we interpolate from the observed
data set by fitting a distribution to those data and then draw from that
distribution for future states in the E variable. Fitting the
distribution is more hypothesis driven because there could be several
different kinds of distribution that could be fitted. Log-normal would
be the classic one and would fit with many datasets and in through use
has the feel of being parsimonius. Log-normal may not be appropriate
though in cases where say there are theoretical true upper limits. Also,
we may wish to capture the possibility of “black swans” that could
dominate what is happening in which case we may wish to use a heavier
tailed distribution. Sean Anderson has done work on this.

Here is a simple example where we use the R library MASS to fit a
distribution to data and then draw future samples from that
distribution. The result of such a sampling after going through your
population model to examine the impacts on a stock will have a smoother
range of possible outcomes.

    library(MASS)
    fittedLN= fitdistr(log(varval$value), densfun="lognormal")
    LN.Esamples= rlnorm(100000, coef(fittedLN)[1], coef(fittedLN)[2])
    plot(density(exp(LN.Esamples)),xlab="E value (not logged)", main="Log-normal distribution fitted")

![](climate.projections_files/figure-markdown_strict/parametric.sampling-1.png)

    # lets try a Weibull distribution which can have a fat tail and may give a higher probability to more extreme events
    fittedWB= fitdistr(varval$value, densfun="weibull")
    WB.Esamples= rweibull(100000, coef(fittedWB)[1], coef(fittedWB)[2])
    plot(density(WB.Esamples),xlab="E value (not logged)", main="Weibull distribution fitted")

![](climate.projections_files/figure-markdown_strict/parametric.sampling-2.png)

    # to my surprise it puts the tail on the lower side which is interesting

What is good about parametric distributions is that you can shift one of
the parameters to immitate a plausible future state. For example, you
can shift the mean up or down or change the variance. Climate change is
projected to often shift the mean and increase the variance and this is
pretty easy to do with most parametric distributions.

Accounting for future climate using down-scaled GCM projections
---------------------------------------------------------------

The IPCC uses &gt;20 global climate models (GCM) to make predictions and
what you often see in IPCC reports in the ensemble median and variance
predictions of this ensemble of models. These down-scaled predictions
are less common for oceanographic variables (they are getting more
common) but they the atmospheric projections are available for many
places. The cliamte atlas of Canada provides a really nice way to obtain
down-scaled predictions for a few carbon emissions scenarios in boxes of
(y\*y) boxes from all over the country. They are not available for areas
that do not touch land though.

What is possible to do is make a statistical model (e.g. a GLM) an
oceanographic variable with the model output during the past (data)
period and then go forward with this. This may not be ideal but at least
it can be better than guessing.

The gslea package has integrated several of the atmospheric climate
projection variables for the Gulf of St. Lawrence already so I will use
that here.

     confint.f= function(x,ylow,yhigh,...){
        polygon(x = c(x, rev(x)), y = c(ylow,rev(yhigh)), border = NA,...)
    }

      proj.var= "Ann.mean.T.med.RCP45"
      tmp= EA.query.f(c(proj.var,"T.deep"),1950:2100,2)
      tmp2= dcast(tmp, year~variable)
      names(tmp2)[3]= "proj.var"
      plot(tmp2$proj.var,tmp2$T.deep)
      rug(tmp2$proj.var)

![](climate.projections_files/figure-markdown_strict/climproject-1.png)

      pred.lm= lm(T.deep~proj.var,data=tmp2)
      summary(pred.lm)

    ## 
    ## Call:
    ## lm(formula = T.deep ~ proj.var, data = tmp2)
    ## 
    ## Residuals:
    ##      Min       1Q   Median       3Q      Max 
    ## -0.86984 -0.24946 -0.09577  0.27351  0.69848 
    ## 
    ## Coefficients:
    ##             Estimate Std. Error t value Pr(>|t|)    
    ## (Intercept)   3.1376     0.5317   5.901 1.16e-06 ***
    ## proj.var      0.5526     0.1238   4.465 8.39e-05 ***
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Residual standard error: 0.3933 on 34 degrees of freedom
    ##   (110 observations deleted due to missingness)
    ## Multiple R-squared:  0.3696, Adjusted R-squared:  0.3511 
    ## F-statistic: 19.94 on 1 and 34 DF,  p-value: 8.387e-05

      tmp2$T.deep.pred= predict(pred.lm,newdata=tmp2)
      plot(tmp2$proj.var,tmp2$T.deep.pred,type="l",col="blue",lwd=3,
           xlab=proj.var, ylab= "Bottom temperature of deep (>200 m) waters EAR 2")
      points(tmp2$proj.var,tmp2$T.deep,pch=20)
      rug(tmp2$proj.var)
      title(main=paste0(proj.var," ensemble median"))

      # so let's create future climate scenarios by calculating the deep water temperature based on the IPCC climate
      # projection and add a residual from the linear model for each to capture the variance in the fitted 
      #relationship
      
      T.deep.pred= predict(pred.lm,newdata=tmp2)+sample(residuals(pred.lm), size=nrow(tmp2), replace=T)
      points(tmp2$proj.var, T.deep.pred, pch=20,col="grey")

![](climate.projections_files/figure-markdown_strict/climproject-2.png)

      #setup a data.frame to output projections
      N=10000
      T.projections= matrix(nrow=nrow(tmp2),ncol=N)
      for (i in 1:N){
        T.projections[,i]= predict(pred.lm,newdata=tmp2)+sample(residuals(pred.lm), size=nrow(tmp2), replace=T)
      }
      
    #In this case the variance carried forward in the future projections is simply the variance in the residuals of the fitted linear model. However one could fit differnt kinds of models any they could also use different GCM projections instead of, or more correctly, in addition to.

    # So lets plot the whole thing with a 90% confidence interval 
    T.proj.quantiles= apply(T.projections,1,quantile,c(0.05, 0.5, 0.95))
    ylims=range(T.proj.quantiles)
    T.proj.quantiles= rbind(tmp2$year,T.proj.quantiles)
        
      plot(T.proj.quantiles[1,],T.proj.quantiles[3,],type="n",ylim=c(ylims[1],ylims[2]),
           xlab="Year", ylab= "Bottom temperature of deep (>200 m) waters")
      confint.f(T.proj.quantiles[1,],T.proj.quantiles[2,],T.proj.quantiles[4,],col="grey")
      lines(T.proj.quantiles[1,],T.proj.quantiles[3,],lwd=3, col="blue")
      title(main=paste0(proj.var," ensemble median"))

![](climate.projections_files/figure-markdown_strict/climproject-3.png)
