#!/usr/bin/env python
# -*- coding: utf-8 -*-

import numpy as np
import csv
from scipy.stats import ttest_rel, ranksums, shapiro
from statsmodels.stats.anova import anova_lm
from scipy import stats
import matplotlib.pyplot as py

f_ind = './RTs_ind_all.csv'
f_group = './RTs_group_all.csv'

def draw_bars(f_ind, f_group):
	compatible_group = []
	compatible_group_err = []
	incompatible_group = []
	incompatible_group_err = []
	with open(f_group, 'rb') as csvfile:
		reader = csv.reader(csvfile, delimiter=',')
		for row in reader:
			if row[0] == 'compatible-go':
				compatible_group.append(float(row[1]))
				compatible_group_err.append(float(row[2]))
			elif row[0] == 'incompatible-go':
				incompatible_group.append(float(row[1]))
				incompatible_group_err.append(float(row[2]))

	compatible_ind = []
	compatible_ind_err = []
	incompatible_ind = []
	incompatible_ind_err = []
	with open(f_ind, 'rb') as csvfile:
		reader = csv.reader(csvfile, delimiter=',')
		for row in reader:
			if row[0] == 'compatible-go':
				compatible_ind.append(float(row[1]))
				compatible_ind_err.append(float(row[2]))
			elif row[0] == 'incompatible-go':
				incompatible_ind.append(float(row[1]))
				incompatible_ind_err.append(float(row[2]))

	g_in = np.mean(incompatible_group_err)*10**(3)/(np.sqrt(len(incompatible_group_err)))
	g_c = np.mean(compatible_group_err)*10**(3)/(np.sqrt(len(compatible_group_err)))
	i_in = np.mean(incompatible_ind_err)*10**(3)/(np.sqrt(len(incompatible_ind_err)))
	i_c = np.mean(compatible_ind_err)*10**(3)/(np.sqrt(len(compatible_ind_err)))
	g_cc = np.mean(compatible_group)*10**(3)
	g_ii = np.mean(incompatible_group)*10**(3)
	i_cc = np.mean(compatible_ind)*10**(3)
	i_ii = np.mean(incompatible_ind)*10**(3)

	py.figure()
	ax = py.subplot(111)
	rects1 = ax.bar([1.3,2.4],[g_ii,i_ii],width=0.4,color='r',yerr=[g_in,i_in],error_kw=dict(elinewidth=1.5,ecolor='black'))
	rects2 = ax.bar([1.7,2.8],[g_cc,i_cc],width=0.4,color='g',yerr=[g_c,i_c],error_kw=dict(elinewidth=1.5,ecolor='black'))
	ax.set_ylabel('Czas reakcji [ms]')
	xTickMarks = ['Grupowy','Indywidualny']
	ax.set_xticks([1.7,2.8])
	xtickNames = ax.set_xticklabels(xTickMarks)
	py.setp(xtickNames, rotation=0, fontsize=20)
	ax.legend((rects1[0], rects2[0]),('Niezgodny', 'Zgodny'))
	ax.set_ylim(340,390)

	py.figure()
	py.boxplot([compatible_group,incompatible_group,compatible_ind,incompatible_ind])
	labels = ('Zgodny-grupowy', 'Niezgodny-grupowy', 'Zgodny-indywidualny', 'Niezgodny-indywidualny')
	py.xticks(range(1,5),labels,rotation=15)
	py.xlabel('Warunki')
	py.ylabel('Czas reakcji [ms]')

	# py.figure()
	# py.hist(compatible_group)
	# py.figure()
	# py.hist(incompatible_group)
	# py.figure()
	# py.hist(compatible_ind)
	# py.figure()
	# py.hist(incompatible_ind)

	py.show()

def test(f_ind, f_group, mode):
	compatible = []
	incompatible = []
	if mode == 'group':
		f_name = f_group
	else:
		f_name = f_ind
	with open(f_name, 'rb') as csvfile:
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
	print 'Median (compatible): ',np.median(compatible)
	print 'Median (incompatible): ',np.median(incompatible)
	### p-value < 0.05 - odrzucamy hipotezę zerową o tym, że brak różnic pomiędzy średnimi
	print "One-way ANOVA p =",stats.f_oneway(compatible,incompatible)[1]  


test(f_ind,f_group,mode='individual')	### mode='group'/'individual'
draw_bars(f_ind,f_group)
