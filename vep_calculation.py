#!/usr/bin/env python
# -*- coding: utf-8 -*-

from obci.analysis.obci_signal_processing import read_manager 
from obci.analysis.obci_signal_processing.tags import tags_file_reader as TFR

# import scipy.stats as stats
import matplotlib.pyplot as py
import numpy as np
from rys_10_20 import rys
import os.path, ast
from scipy.signal import filtfilt, cheby2, butter, resample
import random

def _get_timestamps(tag_file, condition):
	proper_key = _get_proper_key(condition)
	other_condition = 'green' if condition == 'blue' else 'blue'
	keys = ['green-left', 'green-right', 'blue-left', 'blue-right']
	times = {}
	for k in keys:
		times[k] = []
	reader = TFR.TagsFileReader(tag_file)
	tags = reader.get_tags()
	proper_cond = [key for key in keys if condition in key]
	wrong_cond = list(set(keys) - set(proper_cond))
	for tag in tags:
		for c in wrong_cond:
			if tag['name'] == c:
				if isinstance(tag['desc']['key'], basestring) and tag['desc']['key'] != proper_key:
					times[c].append(float(tag['start_timestamp']))
				elif isinstance(tag['desc']['key'], list):
					k = ast.literal_eval(tag['desc']['key'])
					if proper_key not in k:
						times[c].append(float(tag['start_timestamp']))
		for c in proper_cond:
			if tag['name'] == c:
				if isinstance(tag['desc']['key'], basestring) and tag['desc']['key'] == proper_key:
					times[c].append(float(tag['start_timestamp']))
				elif isinstance(tag['desc']['key'], list):
					k = ast.literal_eval(tag['desc']['key'])
					if proper_key in k:
						times[c].append(float(tag['start_timestamp']))				
	return times

def cut_fragments(file_name):
	'''
	format zwracanych danych frags['O1']['blue-right'][[1024],[1024],[1024],...]
	'''
	condition = get_current_condition(file_name)
	mgr = read_manager.ReadManager(file_name+'.xml',
								   file_name+'.raw',
								   file_name+'.tag')
	channels = mgr.get_param('channels_names')
	pointsPerMikroV = mgr.get_param('channels_gains')
	fs = float(mgr.get_param('sampling_frequency'))
	# sig = mgr.get_samples()		### uncomment if signal is already filtered
	sig = _linked_ears_montage(mgr, fs)
	# sig = _one_ear_montage(mgr, fs)
	times = _get_timestamps(file_name+'.tag',condition)
	stim = {}
	for key in times.keys():
		stim[key] = [t*fs for t in times[key]]
	frags_ch = {}
	for ch in channels[:23]:
		frags_ch[ch] = {}
		for key in times.keys():
			frags_ch[ch][key] = []
	for i in xrange(len(sig[:23])):
		sig[i] = sig[i]*float(pointsPerMikroV[i])
		for key in stim.keys():
			for j in stim[key]:
				frags_ch[channels[i]][key].append(sig[i][j:j+int(fs)]-np.mean(sig[i][j-int(0.1*fs):j+int(fs)]))
	frags_ch.pop('A1', None)
	frags_ch.pop('A2', None)
	return frags_ch, fs

def _get_proper_key(condition):
	if condition == 'blue':
		proper_key = 'right'  
	else: 
		proper_key = 'left'
	return proper_key

def get_RTs(file_name):
	condition = get_current_condition(file_name)
	proper_key = _get_proper_key(condition)
	keys = [condition+'-left', condition+'-right']
	rts = {}
	for k in keys:
		rts[k] = []
	tag_file = file_name + '.tag'
	reader = TFR.TagsFileReader(tag_file)
	tags = reader.get_tags()
	for tag in tags:
		if tag['name'] in rts and str(tag['desc']['key']) != 'None':
			if isinstance(tag['desc']['key'], basestring) and tag['desc']['key'] == proper_key:
				rts[tag['name']].append(float(tag['desc']['response_time']))
			else:
				try:
					k = ast.literal_eval(tag['desc']['key'])
					rt = ast.literal_eval(tag['desc']['response_time'])
					ids = [i for i,j in enumerate(k) if j == proper_key]
					rts[tag['name']].append(float(rt[ids[0]]))
				except (TypeError, ValueError) as e:
					print e,' occured in ', file_name
	mean_rts = _compute_mean_RTs(rts)
	return mean_rts

def _filter_signal(signal, fs):
	sig_f = []
	fn = fs/2.
	# b,a = cheby2(1, 10, 0.5/fn, btype = "highpass")
	# b,a = butter(2, np.array([49,51])/fn, btype = "bandstop")
	# d,c = butter(2, 0.1/fn, btype = "highpass")
	d,c = butter(2, np.array([1,12])/fn, btype = "bandpass")
	for i in xrange(len(signal)):
		# sig_f.append(filtfilt(b,a,filtfilt(d, c, signal[i])))
		sig_f.append(filtfilt(d, c, signal[i]))
	return sig_f

def _linked_ears_montage(mgr, fs):
	sig = mgr.get_samples()
	ref_sig = []
	for i in xrange(23):
		new_sig = sig[i]-0.5*(mgr.get_channel_samples('A1')+mgr.get_channel_samples('A2'))
		ref_sig.append(new_sig-np.mean(new_sig))
	filtered_sig = _filter_signal(ref_sig, fs)
	return filtered_sig

def _one_ear_montage(mgr, fs):
	sig = mgr.get_samples()
	ref_sig = []
	for i in xrange(23):
		new_sig = sig[i]-mgr.get_channel_samples('A2')
		ref_sig.append(new_sig-np.mean(new_sig))
	filtered_sig = _filter_signal(ref_sig, fs)
	return filtered_sig

def _compute_mean_RTs(rts):
	mean_rts = {}
	for key in rts.keys():
		data = _reject_outliers(rts[key])
		# print key,'->',list(set(rts[key]) - set(data))	### show differences between lists
		mean_rts[key] = [round(np.mean(data),4), round(np.std(data),4)]
	return mean_rts

def _reject_outliers(data, m=8.):
	d = np.abs(data-np.median(data))
	mdev = np.median(d)
	norm_data = d/mdev
	new_data = [data[i] for i,j in enumerate(norm_data) if j<m]
	return new_data

def get_current_condition(file_name):
	f = os.path.basename(file_name).split('_')
	return 'blue' if f[2].lower() == 'r' or f[2].lower() == 'rb' else 'green'

def plotting(frags, frags_group,condition, ch, mode, fs):
	proper_key = _get_proper_key(condition)
	other_key = 'left' if condition == 'blue' else 'right'

	# if mode == 'no-go':
	# 	condition = 'green' if condition == 'blue' else 'blue'
	'''
	poprawić  dla list!!!!!!!
	'''

	if isinstance(ch, list):
		if mode == 'go':
			rys_go = rys('Go trials')
		else:
			rys_go = rys('No-go trials')
			condition = 'green' if condition == 'blue' else 'blue'
		for ind,i in enumerate(ch):
			av = np.mean(frags[i][condition+'-'+proper_key], 0)
			t = np.linspace(0,len(av)/fs,len(av))
			r1 = rys_go[i].plot(t, av, 'r')
			av_group = np.mean(frags_group[i][condition+'-'+proper_key], 0)
			r2 = rys_go[i].plot(t, av_group, 'g')
			av = np.mean(frags[i][condition+'-'+other_key], 0)
			r3 = rys_go[i].plot(t, av, 'c')
			av_group = np.mean(frags_group[i][condition+'-'+other_key], 0)
			r4 = rys_go[i].plot(t, av_group, 'm')
			for r in [r1, r2]:
				py.setp(r,linewidth=2)
			# rys_go[i].set_xlim(0,0.6)
			# rys_go[i].set_ylim(-6,10)
		rys_go['legend'].plot([],[],'r',label='Individual Compatible')
		rys_go['legend'].plot([],[],'g',label='Group Compatible')
		rys_go['legend'].plot([],[],'c',label='Individual Incompatible')
		rys_go['legend'].plot([],[],'m',label='Group Incompatible')
		rys_go['legend'].legend()
	else:
		fig = py.figure()
		py.title(mode+' trials')
		av = np.array(frags['compatible-'+mode][0])
		av_var = np.array(frags['compatible-'+mode][1])
		t = np.linspace(0,len(av)/fs,len(av))
		av_group = np.array(frags_group['compatible-'+mode][0])
		av_group_var = np.array(frags_group['compatible-'+mode][1])
		py.plot(t, av, c='r', label='Individual Compatible')
		py.fill_between(t, av+av_var, av-av_var, facecolor='r', alpha=0.5)
		py.plot(t, av_group, c='g', label='Group Compatible')
		py.fill_between(t, av_group+av_group_var, av_group-av_group_var, facecolor='g', alpha=0.5)
		av1 = np.array(frags['incompatible-'+mode][0])
		av_var1 = np.array(frags['incompatible-'+mode][1])
		av_group1 = np.array(frags_group['incompatible-'+mode][0])
		av_group_var1 = np.array(frags_group['incompatible-'+mode][1])
		py.plot(t, av1, c='c', label='Individual Incompatible')
		py.fill_between(t, av1+av_var1, av1-av_var1, facecolor='c', alpha=0.5)
		py.plot(t, av_group1, c='m', label='Group Incompatible')
		py.fill_between(t, av_group1+av_group_var1, av_group1-av_group_var1, facecolor='m', alpha=0.5)

		py.legend()
	py.show()

def find_peak(signal, fs):
	deriv = np.gradient(np.array(signal[:0.5*fs]))
	ind = []
	amps = []
	for i in xrange(len(deriv)):
		if i != (0,1):
			if (deriv[i-1] > 0 and deriv[i] < 0):
				ind.append(('max',i-1))
			elif (deriv[i-1] < 0 and deriv[i] > 0):
				ind.append(('min',i-1))
	for i in xrange(len(ind)):
		if i != 0:
			if (ind[i-1][0] == 'min' and ind[i][0] == 'max'):
				amps.append(abs(signal[ind[i][1]] - signal[ind[i-1][1]]))
	return max(amps)

def estimate_mean_variance_bootstrap(frags, channel, fs):
	N = 1000
	amp_cond = {}
	bootstrap_mean = {}
	for c in frags[channel].keys():
		amp = []
		amp_mean = []
		ind = frags[channel][c]
		for i in xrange(N):
			random.shuffle(ind)
			p = int(0.7*len(ind))
			ind_mean = np.mean(ind[:p], axis=0)
			amp_mean.append(ind_mean)
			amp.append(find_peak(ind_mean, fs))
		amp_cond[c] = [np.mean(amp), np.var(amp)]
		bootstrap_mean[c] = [np.mean(amp_mean, axis=0),np.var(amp_mean, axis=0)] 
	return amp_cond, bootstrap_mean


def estimate_mean_variance(frags, channel, fs):
	frags_mean = {}
	for c in frags[channel].keys():
		frags_mean[c] = [np.mean(frags[channel][c], axis=0),np.var(frags[channel][c], axis=0)/len(frags[channel][c])] 
	return frags_mean

def compute_grand_average(f_list_group, f_list_ind):
	keys = ['compatible-go','compatible-no-go','incompatible-go','incompatible-no-go']
	compatible = []
	incompatible = []
	inds = {}
	groups = {}
	for k in keys:
		inds[k] = []
		groups[k] = []
	for i in xrange(len(f_list_group)):
		file_name = './badania_part3/' + f_list_ind[i]
		file_name_group = './badania_part3/' + f_list_group[i]
		times = get_RTs(file_name)
		times1 = get_RTs(file_name_group)

		condition = get_current_condition(file_name)
		proper_key = _get_proper_key(condition)
		other_key = 'left' if condition == 'blue' else 'right'
		other_condition = 'green' if condition == 'blue' else 'blue'

		# compatible.append(times[condition+'-'+proper_key])
		# compatible.append(times1[condition+'-'+proper_key])

		# incompatible.append(times[condition+'-'+other_key])
		# incompatible.append(times1[condition+'-'+other_key])

	# print 'Compatible: ', np.mean(compatible,0)
	# print 'incompatible: ', np.mean(incompatible,0)

		frags, fs = cut_fragments(file_name)
		frags_group, fs = cut_fragments(file_name_group)

		ind = estimate_mean_variance(frags, 'Cz', fs)
		group = estimate_mean_variance(frags_group, 'Cz', fs)
		for j in ind.keys():
			if j == condition+'-'+proper_key:
				inds['compatible-go'].append(ind[j])
				groups['compatible-go'].append(group[j])
			elif j == condition+'-'+other_key:
				inds['incompatible-go'].append(ind[j])
				groups['incompatible-go'].append(group[j])
			elif j == other_condition+'-'+proper_key:
				inds['compatible-no-go'].append(ind[j])
				groups['compatible-no-go'].append(group[j])
			elif j == other_condition+'-'+other_key:
				inds['incompatible-no-go'].append(ind[j])
				groups['incompatible-no-go'].append(group[j])
	f_ind = {}
	f_gr = {}
	for k in inds.keys():
		f_ind[k] = np.mean(inds[k],axis=0)
		f_gr[k] = np.mean(groups[k],axis=0)
	plotting(f_ind, f_gr, 'Cz', 'no-go', fs)
	plotting(f_ind, f_gr, 'Cz', 'go', fs)


if __name__ == "__main__":

	channels = [u'Fp1',u'Fpz',u'Fp2',u'F7',u'F3',u'Fz',u'F4',u'F8',u'T3',u'C3',u'Cz',u'C4',u'T4',u'T5', 
		  u'P3', u'Pz', u'P4', u'T6', u'O1', u'Oz', u'O2']

	path = './badania_part3/'
	f_list_group = ['21_K_RB_01']
	f_list_ind = ['21_K_R_01']

	# compute_grand_average(f_list_group, f_list_ind)

	for i in xrange(len(f_list_ind)):
		condition = get_current_condition(f_list_ind[i])
		frags,fs = cut_fragments(path+f_list_ind[i])
		frags_group,fs = cut_fragments(path+f_list_group[i])
		plotting(frags, frags_group, condition, channels, 'no-go', fs)
		plotting(frags, frags_group, condition, channels, 'go', fs)





