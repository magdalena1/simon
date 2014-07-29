#!/usr/bin/env python
# -*- coding: utf-8 -*-
from __future__ import print_function, division

import numpy as np
from scipy.stats import ttest_rel, ranksums, shapiro
import csv
import glob
from statsmodels.stats.anova import anova_lm
from scipy import stats
from collections import defaultdict
from statsmodels.stats.multicomp import pairwise_tukeyhsd, MultiComparison
import matplotlib.pyplot as py

class ERPStatistics(object):
	def __init__(self,erp_type,channel,condition1=None,condition2=None):

		super(ERPStatistics, self).__init__()
		self.erp_type = erp_type
		f_erp = './csv/poster/'+erp_type+'_'+channel+'_mean.csv'
		self.conditions = self._get_conditions(f_erp)
		self._remove_unknown_values()
		self._test_condition(condition1,condition2)

	def _remove_unknown_values(self):
		for k in self.conditions.keys():
			median = np.median(filter(lambda x: x > 0, self.conditions[k]))
			ids = np.where(np.array(self.conditions[k]) <= 0)
			for i in ids[0]:
				self.conditions[k][i] = median

	def _test_condition(self,condition1,condition2):
		if condition1 in ['group','individual']:
			cond1 = self.conditions['compatible_go_group']+self.conditions['incompatible_go_group']
			cond2 = self.conditions['compatible_go_individual']+self.conditions['incompatible_go_individual']
			self._print_results(cond1,cond2,'group-go','individual-go')
			
			cond1 = self.conditions['compatible_nogo_group']+self.conditions['incompatible_nogo_group']
			cond2 = self.conditions['compatible_nogo_individual']+self.conditions['incompatible_nogo_individual']
			self._print_results(cond1,cond2,'group-no-go','individual-no-go')
		
		elif condition1 in ['compatible','incompatible']:
			cond1 = self.conditions['compatible_nogo_group']+self.conditions['compatible_nogo_individual']
			cond2 = self.conditions['incompatible_nogo_group']+self.conditions['incompatible_nogo_individual']
			self._print_results(cond1,cond2,'compatible-no-go','incompatible-no-go')

			cond1 = self.conditions['compatible_go_group']+self.conditions['compatible_go_individual']
			cond2 = self.conditions['incompatible_go_group']+self.conditions['incompatible_go_individual']
			self._print_results(cond1,cond2,'compatible-go','incompatible-go')

			cond1 = self.conditions['compatible_nogo_group']
			cond2 = self.conditions['incompatible_nogo_group']
			self._print_results(cond1,cond2,'compatible-no-go-group','incompatible-no-go-group')

			cond1 = self.conditions['compatible_go_group']
			cond2 = self.conditions['incompatible_go_group']
			self._print_results(cond1,cond2,'compatible-go-group','incompatible-go-group')

			cond1 = self.conditions['compatible_nogo_individual']
			cond2 = self.conditions['incompatible_nogo_individual']
			self._print_results(cond1,cond2,'compatible-no-go-individual','incompatible-no-go-individual')

			cond1 = self.conditions['compatible_go_individual']
			cond2 = self.conditions['incompatible_go_individual']
			self._print_results(cond1,cond2,'compatible-go-individual','incompatible-go-individual')

	def _print_results(self,cond1,cond2,condition1,condition2):
		print('Mean value (',condition1,'): ',np.mean(cond1))
		print('Mean value (',condition2,'): ',np.mean(cond2))
		print('Shapiro-Wilk test for normality: ',shapiro(cond1)[1],shapiro(cond2)[1]) ### small sample size
		### p-value < 0.05 - mała szansa, że dane pochodzą z rozkładu normalnego
		print('T-test on TWO RELATED samples: ',ttest_rel(cond1,cond2)[1])
		if shapiro(cond1)[1] > 0.05 and shapiro(cond2)[1] > 0.05 and ttest_rel(cond1,cond2)[1] < 0.05:
			print('##################')
		print('Wilcoxon rank-sum statistic: ',ranksums(cond1,cond2)[1])
		if ranksums(cond1,cond2)[1] < 0.05:
			print('!!!!!!!!!!!!!!!!!!')
		### p-value < 0.05 - odrzucamy hipotezę zerową o tym, że brak różnic pomiędzy średnimi - różnice ISTNIEJĄ
		print('One-way ANOVA p =',stats.f_oneway(cond1,cond2)[1])

	def _get_conditions(self,file_name):
		conditions = defaultdict(list)
		with open(file_name, 'rb') as f:
			r = csv.reader(f)
			for row in r:
				if row[0] == 'compatible-go' and row[1] == 'group':
					conditions['compatible_go_group'].append(float(row[2]))
				elif row[0] == 'incompatible-go' and row[1] == 'group':
					conditions['incompatible_go_group'].append(float(row[2]))
				elif row[0] == 'compatible-no-go' and row[1] == 'group':
					conditions['compatible_nogo_group'].append(float(row[2]))
				elif row[0] == 'incompatible-no-go' and row[1] == 'group':
					conditions['incompatible_nogo_group'].append(float(row[2]))
				elif row[0] == 'compatible-go' and row[1] == 'individual':
					conditions['compatible_go_individual'].append(float(row[2]))
				elif row[0] == 'incompatible-go' and row[1] == 'individual':
					conditions['incompatible_go_individual'].append(float(row[2]))
				elif row[0] == 'compatible-no-go' and row[1] == 'individual':
					conditions['compatible_nogo_individual'].append(float(row[2]))
				elif row[0] == 'incompatible-no-go' and row[1] == 'individual':
					conditions['incompatible_nogo_individual'].append(float(row[2]))
		return conditions

if __name__ == '__main__':
	'''
		Possible statistics:
			group vs individual:
				- go
				- no-go
			compatible vs incompatible
				- go
				- no-go
				- go-individual 
				- go-group
				- no-go-individual
				- no-go-group
	'''
	record = np.zeros(1, dtype=[('Electrode','|S12'),('Condition1','|S12'),('go/nogo','|S12'),('Mode','|S12'),('Value','float64')])
	channels = {'prefrontal':['Fp1','Fpz','Fp2'],'frontal':['F3','Fz','F4'],'parietal':['P3','Pz','P4'],'central':['C3','Cz','C4'],'occipital':['O1','Oz','O2']}
	### ('occipital', 'incompatible', 'nogo', 'group', 4.635233271375064)
	# for c in ['Fp1','Fpz','Fp2','C3','C4','Cz','P3','P4','Pz','F3','Fz','F4','O1','Oz','O2']:
	for c in ['Cz','C3','C4']:
		# print('--------------->Channel: ',c)
		s = ERPStatistics('p300',c)
		for k in channels.keys():
			if c in channels[k]:
				electrode = k
		for i,condition in enumerate(s.conditions.keys()):
			header = condition.split('_')
			for value in s.conditions[condition]:
				if header[1] == 'nogo':
					record = np.append(record,np.array([(electrode,header[0],header[1],header[2],value)],dtype=record.dtype))
	record = record[1:]
	# print(record[-1])

	# res2 = pairwise_tukeyhsd(record['Value'],record['Mode'],alpha=0.05)
	# print(res2)
	# mod = MultiComparison(record['Value'],record['go/nogo'])
	# print(mod.groupsunique)
	# results = mod.tukeyhsd()

	g1 = []
	g2 = []
	for i in record:
		if i[3] == 'group':
			g1.append(i[4])
		if i[3] == 'individual':
			g2.append(i[4])

	print(stats.kruskal(g1,g2))
	data = [g1,g2]
	py.boxplot(data,1)
	py.show()

############# nie wychodzą różnice dla wszystkich warunków na raz
	# record_AIO = np.zeros(1, dtype=[('Mode','|S40'),('Value','float64')])
	# for c in ['Fp1','Fpz','Fp2','C3','C4','Cz','P3','P4','Pz','F3','Fz','F4','O1','Oz','O2']:
	# 	# print('--------------->Channel: ',c)
	# 	s = ERPStatistics('p300',c)
	# 	for k in channels.keys():
	# 		if c in channels[k]:
	# 			electrode = k
	# 	for i,condition in enumerate(s.conditions.keys()):
	# 		header = electrode+'_'+condition
	# 		for value in s.conditions[condition]:
	# 			record_AIO = np.append(record_AIO,np.array([(header,value)],dtype=record_AIO.dtype))
	# record_AIO = record_AIO[1:]
	# mod = MultiComparison(record_AIO['Value'],record_AIO['Mode'])
	# results = mod.tukeyhsd()
#############

	# s = ERPStatistics('p300','Cz','group','individual')

	# res2 = pairwise_tukeyhsd(record['Value'],record['go/nogo'])
	# print(res2)
	# py.plot([0,1,2], res2[1][2], 'o')
	# py.errorbar([0,1,2], res2[1][2], yerr=np.abs(res2[1][4].T-res2[1][2]), ls='o')
	# xlim = -0.5, 2.5
	# py.hlines(0, *xlim)
	# py.xlim(*xlim)
	# pair_labels = mod.groupsunique[np.column_stack(res2[1][0])]
	# py.xticks([0,1,2], pair_labels)
	# py.title('Multiple Comparison of Means - Tukey HSD, FWER=0.05' + '\n Pairwise Mean Differences')
	# py.show()