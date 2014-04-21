#!/usr/bin/env python
# -*- coding: utf-8 -*-

import numpy as np
from scipy.stats import ttest_rel, ranksums, shapiro
import csv

f_erp = './p300_Oz.csv'
condition_1 = 'compatible-go'
condition_2 = 'incompatible-go'

condition1 = []
condition2 = []
with open(f_erp, 'rb') as f:
	r = csv.reader(f)
	for row in r:
		if row[0] == condition_1:
			condition1.append(float(row[2]))
		elif row[0] == condition_2:
			condition2.append(float(row[2]))

print 'Shapiro-Wilk test for normality: ',shapiro(condition1)[1],shapiro(condition2)[1] ### small sample size
### jeżeli p-value < 0.05 - mała szansa, że dane pochodzą z rozkładu normalnego
print 'Mean value (',condition_1,'): ',np.mean(condition1)
print 'Mean value (',condition_2,'): ',np.mean(condition2)
print 'T-test on TWO RELATED samples: ',ttest_rel(condition1,condition2)[1]
print 'Wilcoxon rank-sum statistic: ',ranksums(condition1,condition2)[1]
### p-value < 0.05 - odrzucamy hipotezę zerową o tym, że brak różnic pomiędzy średnimi - różnice ISTNIEJĄ

