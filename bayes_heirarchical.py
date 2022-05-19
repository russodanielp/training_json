
"""
Build a heirarchical bayesian model

Here we to estimate the positive predictive and sensitivity of different gene targets available
for pubchem.

This seems to be a good resource....
https://docs.pymc.io/projects/examples/en/latest/case_studies/hierarchical_partial_pooling.html

Heres another maybe better resource:
https://blog.dominodatalab.com/ab-testing-with-hierarchical-models-in-python

"""


import pandas as pd
import pymc3 as pm
import numpy as np



def do_bayes(data: pd.DataFrame):
    """ builds a heirarchical model for every gene in the data.  Columns should be
     GeneSymbol, TP, FP, FN, and TN counts """


    ## using a uniform prior for the probabilities
    # was converging around 0.5.  Not sure why.
    # following this: https://blog.dominodatalab.com/ab-testing-with-hierarchical-models-in-python
    # it suggests using a beta distribution

    with pm.Model() as heir_model:

        a = pm.Exponential("a", 1)
        b = pm.Exponential("b", 1)

        # priors
        p_sens = pm.Beta("p_sens", a, b, shape=data.shape[0])
        p_ppv = pm.Beta("p_ppv", a, b, shape=data.shape[0])

        # Set of observations, in this case we have two observation datasets.
        # liklihoods
        sens = pm.Binomial("obs_sens", p=p_sens, n=data.TP+data.FN, observed=data.TP)
        ppv = pm.Binomial("obs_ppv", p=p_ppv, n=data.TP+data.FP, observed=data.TP)


    with heir_model:
        step = pm.Metropolis()
        trace = pm.sample(draws=1000, tune=2000, target_accept=0.95, return_inferencedata=False, cores=1)


    pd.DataFrame(trace["p_sens"], columns=data.GeneSymbol).to_csv('data/p_sens.csv')
    pd.DataFrame(trace["p_ppv"], columns=data.GeneSymbol).to_csv('data/p_ppv.csv')

if __name__ == '__main__':

    data = pd.read_csv('data/for_bayes.csv')
    print(data)
    do_bayes(data)