## PAML site-model comparison

In this analysis the 2 nested site models - NS7 (also called M7; estimates beta paramters) and NS8 (also called M8; estimates beta & w>1 parameters) are compared.
These scripts are run using command `codeml NS7/codonml_site_model_NS7.ctl` and `codeml NS8/codonml_site_model_NS8.ctl`

The fits between the models are compared by plugging their maximum likelihoods into a Likelihood Ratio Test: LRT= lnL(M8)-lnL(M7). 
Compare 2*LRT to Χ2 to calculate p-value (df=2).
