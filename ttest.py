#!/usr/bin/env python
# -*- coding: utf-8 -*-

import numpy as np
from scipy.stats import ttest_rel, ranksums, shapiro

individual_incompatible = [0.8099,3.5456,10.0342,3.7482,7.3838,4.5038,20.8913,10.8632,8.9243]
group_incompatible = [4.8932,8.9152,12.0748,8.3127,9.4891,7.6284,16.0455,11.5416,7.8138]
individual_compatible = [0.9407,7.4532,9.8394,2.698,5.6294,8.5416,19.3103,10.9919,5.5789]
group_compatible = [5.0762,12.3490,13.6093,10.3565,6.7808,11.1409,21.3137,10.5191,11.7236]

individual = individual_compatible+individual_incompatible
group = group_compatible+group_incompatible

print 'Shapiro-Wilk test for normality: ',shapiro(individual)[1],shapiro(group)[1] ### small sample size
### jeżeli p-value < 0.05 - mała szansa, że dane pochodzą z rozkładu normalnego
print 'Mean value (individual): ',np.mean(individual)
print 'Mean value (group): ',np.mean(group)
print 'T-test on TWO RELATED samples: ',ttest_rel(individual,group)[1]
print 'Wilcoxon rank-sum statistic: ',ranksums(individual,group)[1]
### p-value < 0.05 - odrzucamy hipotezę zerową o tym, że brak różnic pomiędzy średnimi - różnice ISTNIEJĄ
