#!/usr/bin/env python
# -*- coding: utf-8 -*-

''' generuje 4 rodzaje tagów '''

from obci.analysis.obci_signal_processing.tags import tags_file_writer as tags_writer
from obci.analysis.obci_signal_processing.tags import tag_utils
from obci.analysis.obci_signal_processing import read_manager 
import csv, ast
import matplotlib.pyplot as py

def find_blinking_dioda(dioda,fs):
	a, n = 0, 0
	blinks = []
	syg_dioda = dioda
	# py.plot(syg_dioda)
	while n < len(syg_dioda):
		if syg_dioda[n] < 0.8:
			if a == 0:
				a = 1
		if syg_dioda[n] > 0.95:
			if a == 1:
				a = 0
				blinks.append(n)
				py.axvline(n,c='r')
		n += 1
	
	# py.show()
	return blinks

def read_csv(csvfile):
	rows = []
	with open(csvfile,'rb') as f:
		csv_reader = csv.reader(f, delimiter=',')
		header = csv_reader.next()
		key_id = [i for i,j in enumerate(header) if j == 'trial_key.keys']
		rt_id = [i for i,j in enumerate(header) if j == 'trial_key.rt']
		cond_id = [i for i,j in enumerate(header) if j == 'condition']
		pos_id = [i for i,j in enumerate(header) if j == 'stimulus_position']
		for row in csv_reader:
			if row[0]:
				params = {}
				params['color'] = row[0].split('.')[0]
				params['position'] = 'left' if float(row[pos_id[0]]) < 0 else 'right'
				if row[key_id[0]] != 'None':
					keys = ast.literal_eval(row[key_id[0]])
					if len(keys) > 1:
						params['key'] = []
						for item in keys:
							params['key'].append('right') if item == 'p' else params['key'].append('left')
					else:
						params['key'] = 'right' if keys[0] == 'p' else 'left'
				else:
					params['key'] = None
				if row[rt_id[0]]:
					rts = ast.literal_eval(row[rt_id[0]])
					if len(rts) > 1:
						params['response_time'] = []	
						for item in rts:
							params['response_time'].append(item)
					else:
						params['response_time'] = rts[0]
				else:
					params['response_time'] = ''
				params['condition'] = row[cond_id[0]]
				rows.append(params)
	return rows

def tag_writer(fs, dioda, csvfile, file_name):
	csv_params = read_csv(csvfile)
	blinks = find_blinking_dioda(dioda/dioda.max(),fs)
	writer = tags_writer.TagsFileWriter(file_name+'.tag')
	# print file_name,len(csv_params),len(blinks)
	for i in xrange(len(csv_params)):
		tag = tag_utils.pack_tag_to_dict(blinks[i]/fs, blinks[i]/fs+0.75, 
										 csv_params[i]['color']+'-'+csv_params[i]['position'], {
										'response_time': csv_params[i]['response_time'],
										'key': csv_params[i]['key'],
										'condition': csv_params[i]['condition']})
		writer.tag_received(tag)
	writer.finish_saving(0.0)

if __name__ == "__main__":
	f_list = ['36_K_LB_01','36_K_L_01']
	csv_list = ['35_K_RB_01_2014_Mar_12_1214','36_K_L_01_2014_Mar_12_1226']
	path = './badania_part5/'

	for i,f in enumerate(f_list):
		file_name = path+f
		csvfile = path+csv_list[i]+'.csv'
		mgr = read_manager.ReadManager(file_name+'.xml',
								   		file_name+'.raw',
								   		file_name+'.tag')

		dioda = mgr.get_channel_samples('diode')
		fs = float(mgr.get_param('sampling_frequency'))
		tag_writer(fs, dioda, csvfile, file_name)
