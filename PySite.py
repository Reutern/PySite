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

baseDic = {'A':0, 'a':0, 'C':1, 'c':1, 'G':2, 'g':2, 'T':3, 't':3}

# A simple function that reads the input sequences in the FASTA format
def read_fasta(seqFile):
	seqs = []
	seqNames = []
	f_in = open(seqFile, 'r')
	while(True):
		# Read name
		line_name = f_in.readline()
		if not line_name:
			break 		# End of file
		elif line_name[0] != '>':
			continue	# Line does not belong to a sequence. Read new line!
		elif line_name[0] == '>':
			seqNames.append(line_name.strip()[1:])

		# Read sequence 
		line_seq = f_in.readline().strip()

		# Check sequence 
		for base in line_seq:
			if(base not in ['A', 'a', 'C', 'c', 'G', 'g', 'T', 't']):
				print("Error in file {0}: {1} is not a legal base!".format(seqFile, base))
				return ValueError
		seqs.append(line_seq)

	assert(len(seqNames) == len(seqs))
	f_in.close()
	return (seqs, seqNames)


# A function that reads the motif file
def read_wtmx(motifFile, pseudo_counts):
	motifs = []
	motifNames = []
	f_in = open(motifFile, 'r')
	while(True):
		# Read name
		line_name = f_in.readline()
		if not line_name:
			break 		# End of file
		elif line_name[0] != '>':
			continue		# Line does not belong to a sequence. Read new line!
		elif line_name[0] == '>':
			(name_tmp, length) = line_name[1:].split()
			motifNames.append(name_tmp)
				
		pwm = np.zeros([int(length), 4])
		# Read pwm positions
		for pos in range(int(length)):
			line_pwm = f_in.readline().split()
			if len(line_pwm) != 4:
				print("Error in file {0}: motif {1} has corrupt position {2}".format(motifFile, name_tmp, pos))
				return ValueError
			
			for idx in range(4):
				pwm[pos, idx] = float(line_pwm[idx])
		# Apply pseudo counts
		if pseudo_counts > 0:
			sum_row = sum(pwm[0])
			pwm += sum_row * pseudo_counts
		motifs.append(pwm)
	assert(len(motifNames) == len(motifs))
	f_in.close()

	return (motifs, motifNames)


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
		else:
			print("Error in function: reverse_complement! {0} is not a legal base".format(base))
			return ValueError
		seq_reverse += base_reverse
	assert(len(seq) == len(seq_reverse))
	return seq_reverse


# A function for reading the sites of a single motif in a sequnce
def read_sites(seq, motif):
	seqLength = len(seq)
	motifLength = len(motif)
	Npos = seqLength - motifLength + 1

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
		seq_tmp = seq[pos:pos+motifLength]
		seq_tmp_reverse = reverse_complement(seq_tmp)
		for idx in range(motifLength):
			base_plus = baseDic[seq_tmp[idx]]
			weight_plus *= motif[idx, base_plus]
			base_reverse = baseDic[seq_tmp_reverse[idx]]
			weight_reverse *= motif[idx, base_reverse]
		sites_plus[pos] = weight_plus / weight_consensus		
		sites_reverse[pos] = weight_reverse / weight_consensus		
	return (sites_plus, sites_reverse)
			
 
def main(argv=None):

	seqFile = ''
	motifFile = ''
	outputDir = '.'
	backFile = ''
	pseudo_counts = 0

	try:
		opts, args = getopt.getopt(argv,"hs:m:c:o:b:",[])
	except getopt.GetoptError:
		print 'PySite.py -s <input_seq> -m <input_motif> -o <output_dir> [-c <pseudo_counts def=0.0> -b <background_track>]'
		sys.exit(2)

	for opt, arg in opts:
		if opt == '-h':
			print 'PySite.py -s <input_seq> -m <input_motif> -o <output_dir> [-c <pseudo_counts def=0.0> -b <background_track>]'
			sys.exit()
		elif opt in ("-s"):
			seqFile = arg
		elif opt in ("-m"):
			motifFile = arg
		elif opt in ("-c"):
			pseudo_counts = float(arg)		
		elif opt in ("-o"):
			outputDir = arg
		elif opt in ("-b"):
			backFile = arg
		

	# Test whether input files exist
	if not os.path.isfile(seqFile):
		sys.exit(seqFile + " is not a file!")
	if not os.path.isfile(motifFile):
		sys.exit(motifFile + " is not a file!")

	# Read background track File if existing
	background = ''
	if os.path.isfile(backFile):
		background = np.genfromtxt(backFile)	
	try: 
		(seqs, seqNames) = read_fasta(seqFile)
		(motifs, motifNames) = read_wtmx(motifFile, pseudo_counts)
	except ValueError:
		print "Reading Error"
		sys.exit()
	scores_avg = np.zeros(142)
	for idx_seq, seq in enumerate(seqs):
		for idx_motif, motif in enumerate(motifs):
			(sites_plus, sites_reverse) = read_sites(seq, motif)
			Nsites = len(sites_plus)
			range_idx = range(1, Nsites + 1)
			fig, ax = pylab.subplots(ncols=1)
			if background != '':
				ax2 = ax.twinx()
				length = background[idx_seq].size
				range_bg = range(1, length + 1)
				ax2.plot(range_bg, background[idx_seq], linestyle='-', linewidth=2, color = 'k', alpha = 0.5, 
						label="MNase dyad score max = " + str("%.3g" % max(background[idx_seq])))		
			ax.plot(range_idx, sites_reverse, 
                linestyle='-', linewidth=4, alpha = 0.75, label="- strand  max = " + str("%.3g" % max(sites_reverse)))	
			ax.plot(range_idx, sites_plus, 
                linestyle='-', linewidth=4, alpha = 0.75, label="+ strand  max = " + str("%.3g" % max(sites_plus)))
			
			# Show the strongest binding sites
			fig.gca().set_position((.1, .3, .8, .6))
			motif_length = len(motif)			
			best_plus = np.argmax(sites_plus)
			seq_flank_1 = seq[max(best_plus-10, 0):best_plus]
			seq_site = seq[best_plus:best_plus+motif_length]
			seq_flank_2 = seq[best_plus+motif_length:min(best_plus+motif_length+10, len(seq))]
			fig.text(.1, .1, '+ strand:  ' + seq_flank_1 + '_' + seq_site + '_' + seq_flank_2)

			best_reverse = np.argmax(sites_reverse)
			seq_flank_1 = seq[max(best_reverse-10, 0):best_reverse]
			seq_site = seq[best_reverse:best_reverse+motif_length]
			seq_flank_2 = seq[best_reverse+motif_length:min(best_reverse+motif_length+10, len(seq))]
			fig.text(.1, .05, ' -  strand:  ' + seq_flank_1 + '_' + seq_site + '_' + seq_flank_2)


			title_tmp = "{0} sites in {1}".format(motifNames[idx_motif], seqNames[idx_seq])
			plt.title(title_tmp, size = 18)
   			ax.legend(loc='upper right', fontsize = 10)
			if background != '':
	   			ax2.legend(loc='upper left', fontsize = 10)
			ax.set_xlabel("position")
			ax.set_ylabel("binding weight")
			if background != '':
				ax2.set_ylabel("MNase dyad score")
			fig_name = outputDir + '/sites_' + seqNames[idx_seq] + '_' + motifNames[idx_motif] + '.png' 
			fig.savefig(fig_name, bbox_inches='tight')


if __name__ == '__main__':
	main(sys.argv[1:])  




