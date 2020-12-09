README by Zane Kliesmete (21.04.2020)

In this analysis the 2 nested site models - NS7 (also called M7; estimates beta paramters) and NS8 (also called M8; estimates beta & w>1 parameters) are compared.
Each is run using command 'codeml script.ctl'

The fits between the models are compared by plugging their maximum likelihoods in a Likelihood Ratio Test: LRT= lnL(M8)-lnL(M7). 
Compare 2*LRT to Χ2 to calculate p-value (df=2).