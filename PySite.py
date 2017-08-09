#!/usr/bin/env python
#
#Created 13.10.2016
#
#@author: Marc von Reutern
#
#
#Programm for the detection of binding sites

import pylab
import matplotlib.pyplot as plt
import numpy as np
import sys, getopt
import os.path
plt.style.use('ggplot')
plt.rcParams['axes.facecolor']='w'

baseDic = {'A':0, 'a':0, 'C':1, 'c':1, 'G':2, 'g':2, 'T':3, 't':3, 'N':10, 'n':10, '-':100}

# A simple function that reads the input sequences in the FASTA format
def read_fasta(seq_file):
	seqs = []
	seq_names = []
	f_in = open(seq_file, 'r')
	while(True):
		# Read name
		line_name = f_in.readline()
		if not line_name:
			break 		# End of file
		elif line_name[0] != '>':
			continue	# Line does not belong to a sequence. Read new line!
		elif line_name[0] == '>':
			seq_names.append(line_name.strip()[1:])

		# Read sequence 
		line_seq = f_in.readline().strip()

		# Check sequence 
		for base in line_seq:
			if(base not in baseDic.keys()):
				print("Error in file {0}: {1} is not a legal base!".format(seq_file, base))
				return ValueError
		seqs.append(line_seq)

	assert(len(seq_names) == len(seqs))
	f_in.close()
	return (seqs, seq_names)


# A function that reads the motif file
def read_wtmx(motif_file, pseudo_counts):
	motifs = []
	motif_names = []
	f_in = open(motif_file, 'r')
	while(True):
		# Read name
		line_name = f_in.readline()
		if not line_name:
			break 		# End of file
		elif line_name[0] != '>':
			continue		# Line does not belong to a sequence. Read new line!
		elif line_name[0] == '>':
			(name_tmp, length) = line_name[1:].split()
			motif_names.append(name_tmp)
				
		pwm = np.zeros([int(length), 4])
		# Read pwm positions
		for pos in range(int(length)):
			line_pwm = f_in.readline().split()
			if len(line_pwm) != 4:
				print("Error in file {0}: motif {1} has corrupt position {2}".format(motif_file, name_tmp, pos))
				return ValueError
			
			for idx in range(4):
				pwm[pos, idx] = float(line_pwm[idx])
		# Apply pseudo counts
		if pseudo_counts > 0:
			sum_row = sum(pwm[0])
			pwm += sum_row * pseudo_counts
		motifs.append(pwm)
	assert(len(motif_names) == len(motifs))
	f_in.close()

	return (motifs, motif_names)


# A function to build the reverse complement of a single sequence
def reverse_complement(seq):
	seq_reverse = ''
	for pos in range(len(seq)):
		base = seq[-pos-1]
		base_reverse = ''
		if base in ['A', 'a']:
			base_reverse = 'T'
		elif base in ['C', 'c']:
			base_reverse = 'G'
		elif base in ['G', 'g']:
			base_reverse = 'C'
		elif base in ['T', 't']:
			base_reverse = 'A'
		elif base in ['N', 'n']:
			base_reverse = 'N'
		elif base == '-':
			base_reverse = '-'
		else:
			print("Error in function: reverse_complement! {0} is not a legal base".format(base))
			return ValueError
		seq_reverse += base_reverse
	assert(len(seq) == len(seq_reverse))
	return seq_reverse


# A function for reading the sites of a single motif in a sequnce
def read_sites(seq, motif, padding=True):

	# Check the motif input. If motif is a list, check if the entry looks like a motif
	if type(motif) is list:
		if len(motif) == 1 and type(motif[0]) is np.ndarray and motif[0].shape[1] == 4:		
			motif = motif[0]
		else:
			sys.exit("read_sites() takes np.ndarray as motif input, instead of a list with {} entries!".format(len(motif)))			

	if motif.shape[1] != 4:
		sys.exit("read_sites(): motif input corrupt. Too many columns: {}".format(motif.shape[1]))					

	seqLength = len(seq)
	motifLength = len(motif)
	Npos = seqLength + motifLength + 1

	# Extend the sequence on both sites with 'N' or '-' if padding is disabled
	if padding:
		seq_padded = 'N' * motifLength + seq + 'N' * motifLength
	else:
		seq_padded = '-' * motifLength + seq + '-' * motifLength

	# get the reference score for the consensus site
	weight_consensus = 1
	for idx in range(motifLength):
		weight_consensus *= max(motif[idx,:])

	# go through all positions of the plus and reverse strand
	sites_plus = np.zeros(Npos)
	sites_reverse = np.zeros(Npos)
	for pos in range(Npos):
		weight_plus = 1
		weight_reverse = 1
		seq_tmp = seq_padded[pos:pos+motifLength]
		seq_tmp_reverse = reverse_complement(seq_tmp)
		for idx in range(motifLength):
			base_plus = baseDic[seq_tmp[idx]]
			if base_plus < 4:
				weight_plus *= motif[idx, base_plus]
			elif base_plus == 10:
				weight_plus *= min(motif[idx,:])	# if base is a N, score worst possible option
			elif base_plus == 100:
				weight_plus = 0						# if position is missing (base = '-'), score no binding
			base_reverse = baseDic[seq_tmp_reverse[idx]]
			if base_reverse < 4:
				weight_reverse *= motif[idx, base_reverse]
			elif base_reverse == 10:
				weight_reverse *= min(motif[idx,:])	# if base is a N, score worst possible option
			elif base_reverse == 100:
				weight_reverse = 0					# if position is missing (base = '-'), score no binding
		sites_plus[pos] = weight_plus / weight_consensus		
		sites_reverse[pos] = weight_reverse / weight_consensus		
		
	return (sites_plus, sites_reverse)
			
 
# Parameter reading function (from file)
def read_parameters(parameter_file):
	params = {} 
	if not os.path.isfile(parameter_file):
		sys.exit(parameter_file + " is not a file!")
	f_in = open(parameter_file, 'r')
	for line in f_in.readlines():
		if '=' not in line:	# comment or empty line
			continue
		line = line.split('=')
	   	key_tmp = line[0].strip()
		value_tmp = ''
		if len(line) > 2:
	   		sys.exit("Wrong format! Can not read " + parameter_file)	 
		if not line[1].isspace():
		   	value_tmp = line[1].split()[0].strip()	# second split to get rid of comments	   			
		params.update({key_tmp : value_tmp})
	f_in.close()
	return params 


# Parameter writing function (to file)
def save_parameters(parameter_file, params):
	f_in = open(parameter_file, 'w')
	for par in params.keys():
		f_in.write(par + ' = ' + str(params[par]) + '\n')
	f_in.close()
	return params 

 
def main(argv=None):
	seq_file = ''
	motif_file = ''
	output_dir = '.'
	site_threshold = 0.05
	background_file = ''
	pseudo_counts = 0
	padding = True

	# Parameter file and reading conditions
	parameter_file = 'param.txt'
	read_parameter_file = False
	save_parameter_file = True

	doc_string = 'PySite.py -s <input_seq> -m <input_motif> -o <output_dir> '
	doc_string_optional = ' [-t <site_threshold> -p <parameter_file> -c <pseudo_counts def=0.0> -b <background_track> -u (disable seq. padding)]'

	try:
		opts, args = getopt.getopt(argv,"hs:m:c:o:t:b:p:u",[])
	except getopt.GetoptError:
		print doc_string + doc_string_optional
		sys.exit(2)
	for opt, arg in opts:
		if opt == '-h':
			print doc_string + doc_string_optional
			sys.exit()
		elif opt in ("-p"):
			parameter_file = arg
			read_parameter_file = True
			save_parameter_file = False
			break
		elif opt in ("-s"):
			seq_file = arg
		elif opt in ("-m"):
			motif_file = arg
		elif opt in ("-c"):
			pseudo_counts = float(arg)		
		elif opt in ("-o"):
			output_dir = arg
		elif opt in ("-t"):
			site_threshold = float(arg)		
		elif opt in ("-b"):
			background_file = arg
		elif opt in ("-u"):
			padding = False

	# Check if all important parameters are set
	if(seq_file == '' or motif_file == '' or output_dir == ''):
		print 'Reading Parameters'
		read_parameter_file = True
		save_parameter_file = False

	# Read parameter file	
	if (read_parameter_file):
		params = read_parameters(parameter_file) 
		try: 
			seq_file = params['seq_file']
			motif_file = params['motif_file']
			output_dir = params['output_dir']
		except KeyError:
			sys.exit("Reading Error")
		try:
			background_file = params['background_file']
		except KeyError:
			pass	# do not update predefined parameter
		try:
			padding = params['padding']
		except KeyError:
			pass	# do not update predefined parameter
		try:
			site_threshold = float(params['site_threshold'])
		except (KeyError, ValueError):
			pass	# do not update predefined parameter
		try:
			pseudo_counts = float(params['pseudo_counts'])
		except (KeyError, ValueError):
			pass	# do not update predefined parameter
			
	# Save parameters if not specified otherwise
	if (save_parameter_file):
		params = {}
	   	params['seq_file'] = seq_file
		params['motif_file'] = motif_file
		params['output_dir'] = output_dir 
		params['background_file'] = background_file 
		params['site_threshold'] = site_threshold 		
		params['pseudo_counts'] = pseudo_counts 
		params['padding'] = padding
		params = save_parameters(parameter_file, params) 

	# Test whether input files exist
	if not os.path.isfile(seq_file):
		sys.exit(seq_file + " is not a file!")
	if not os.path.isfile(motif_file):
		sys.exit(motif_file + " is not a file!")

	# Read background track file if existing
	background = ''
	if os.path.isfile(background_file):
		background = np.genfromtxt(background_file)	
	try: 
		(seqs, seq_names) = read_fasta(seq_file)
		(motifs, motif_names) = read_wtmx(motif_file, pseudo_counts)
	except ValueError:
		print "Reading Error"
		sys.exit()

	# go through all sequences and motifs	
	for idx_seq, seq in enumerate(seqs):
		for idx_motif, motif in enumerate(motifs):
			(sites_plus, sites_reverse) = read_sites(seq, motif, padding)
			Nsites = len(sites_plus)
			motif_len = len(motif)
			range_idx = range(-motif_len+1, Nsites + 1 - motif_len)

			# write list of binding sites to a file
			f_so = open(output_dir + '/sites_' + seq_names[idx_seq] + '_' + motif_names[idx_motif] + '.txt' ,'w')
			if padding:
				seq_padded = ('N' * motif_len + seq + 'N' * motif_len)
			else:
				seq_padded = ('-' * motif_len + seq + '-' * motif_len)
			for idx_pos in range(1, len(sites_plus)-1):
				seq_tmp = seq_padded[idx_pos:idx_pos+motif_len]
				seq_rev_tmp = reverse_complement(seq_tmp)
				if sites_plus[idx_pos] > site_threshold:
					f_so.write(str(idx_pos-motif_len+1) + '\t' + seq_tmp + '\t+\t' + str(sites_plus[idx_pos]) + '\n')
				if sites_reverse[idx_pos] > site_threshold:					
					f_so.write(str(idx_pos-motif_len+1) + '\t' + seq_rev_tmp + '\t-\t' + str(sites_reverse[idx_pos]) + '\n')
			f_so.close()								
								
			# plot background and bidning site composition
			fig, ax = pylab.subplots(ncols=1)
			if background != '':
				ax2 = ax.twinx()
				length = background[idx_seq].size
				range_bg = range(1, length + 1 - 2 * len(motif))
				ax2.plot(range_bg, background[idx_seq], linestyle='-', linewidth=2, color = 'k', alpha = 0.5, 
						label="background = " + str("%.3g" % max(background[idx_seq])))		
			ax.plot(range_idx, sites_reverse, 
				linestyle='-', linewidth=4, alpha = 0.75, label="- strand  max = " + str("%.3g" % max(sites_reverse)))	
			ax.plot(range_idx, sites_plus, 
				linestyle='-', linewidth=4, alpha = 0.75, label="+ strand  max = " + str("%.3g" % max(sites_plus)))
			
			# show the strongest binding sites
			fig.gca().set_position((.1, .3, .8, .6))
			motif_length = len(motif)			
			best_plus = np.argmax(sites_plus) - motif_length
			seq_flank_1 = seq[max(best_plus-10, 0):best_plus]
			seq_site = seq[best_plus:best_plus+motif_length]
			seq_flank_2 = seq[best_plus+motif_length:min(best_plus+motif_length+10, len(seq))]
			fig.text(.1, .1, '+ strand:  ' + seq_flank_1 + '_' + seq_site + '_' + seq_flank_2)

			best_reverse = np.argmax(sites_reverse) - motif_length
			seq_flank_1 = seq[max(best_reverse-10, 0):best_reverse]
			seq_site = seq[best_reverse:best_reverse+motif_length]
			seq_flank_2 = seq[best_reverse+motif_length:min(best_reverse+motif_length+10, len(seq))]
			fig.text(.1, .05, ' -  strand:  ' + seq_flank_1 + '_' + seq_site + '_' + seq_flank_2)

			title_tmp = "{0} sites in {1}".format(motif_names[idx_motif], seq_names[idx_seq])
			plt.title(title_tmp, size = 18)
   			ax.legend(loc='upper right', fontsize = 10)
   			ax.set_xlim([-len(motif)+1, Nsites + 1 - len(motif)])
			if background != '':
	   			ax2.legend(loc='upper left', fontsize = 10)
			ax.set_xlabel("position")
			ax.set_ylabel("binding weight")
			if background != '':
				ax2.set_ylabel("background")
			fig_name = output_dir + '/sites_' + seq_names[idx_seq] + '_' + motif_names[idx_motif] + '.png' 
			fig.savefig(fig_name, bbox_inches='tight')


if __name__ == '__main__':
	main(sys.argv[1:])  





