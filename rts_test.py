#!/usr/bin/env python
# -*- coding: utf-8 -*-

import numpy as np
import csv
from scipy.stats import ttest_rel, ranksums, shapiro

f_ind = './RTs_ind_all.csv'
f_group = './RTs_group_all.csv'

compatible = []
incompatible = []
for f in [f_ind, f_group]:
	with open(f, 'rb') as csvfile:
		reader = csv.reader(csvfile, delimiter=',')
		for row in reader:
			if row[0] == 'compatible-go':
				compatible.append(float(row[1]))
			elif row[0] == 'incompatible-go':
				incompatible.append(float(row[1]))
print 'T-test on TWO RELATED samples: ',ttest_rel(compatible,incompatible)[1]
print 'Wilcoxon rank-sum statistic: ',ranksums(compatible,incompatible)[1]
print 'Shapiro-Wilk test for normality: ',shapiro(compatible)[1],shapiro(incompatible)[1] ### small sample size
### jeżeli p-value < 0.05 - mała szansa, że dane pochodzą z rokładu normalnego
print 'Mean value (compatible): ',np.mean(compatible)
print 'Mean value (incompatible): ',np.mean(incompatible)
### p-value < 0.05 - odrzucamy hipotezę zerową o tym, że brak różnic pomiędzy średnimi
