

from scitbx.math import gamma_complete, gamma_incomplete, gamma_incomplete_complement
from scitbx.math import distributions

def chi_sq_test(vals):
    '''Does a chi-squared test on sum(xi**2) to find the probability of pulling a chi-squared result less extreme than this one'''

    chi_sq = sum([x**2 for x in vals])
    k = len(vals)

    return 1-gamma_incomplete(k/2, chi_sq/2)

def convert_pvalue_to_zscore(pval, two_tailed=True):
    """Convert p-value back to a z-score"""

    # If two-tailed test, need to halve the p-value
    if two_tailed:
        pval = pval/2

    # Prob of being not in pval region
    prob = 1-pval

    # Create normal distribution to convert - N(0,1)
    nrm = distributions.normal_distribution()

    return nrm.quantile(prob)


vals = [1.1]*100
print vals[0]
pval = chi_sq_test(vals)
print pval
print convert_pvalue_to_zscore(pval)

vals = [1.4]*100
print vals[0]
pval = chi_sq_test(vals)
print pval
print convert_pvalue_to_zscore(pval)




