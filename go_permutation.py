#!/usr/bin/python

from scipy import stats as ss
import numpy
import statsmodels.stats.multitest  as multitest
import random
import time
import sys

thread_num = int(sys.argv[1])
num_permut = 0

if thread_num <= 4:
	num_permut = 1429
else:
	num_permut = 1428

in_file = open("../promoter_go_terms/go_term_abundance_by_gene.table","r")

go_terms = []
gene_list = []
abundance = []
first = True
for line in in_file:
	line = line.strip("\n")
	frags = line.split("\t")
	if first:
		first = False
		go_terms = frags[1:]
	else:
		gene_list.append(frags[0])
		abundance.append(frags[1:])

in_file.close()
abundance = numpy.asarray(abundance,dtype = int)


in_file = open("../promoter_go_terms/total_tf_by_gene.table","r")
tf_total = []
first = True
for line in in_file:
	if first:
		first = False
	else:
		tf_total.append(int(line.strip("\n").split("\t")[1]))
temp_list = range(5)

in_file.close()
out_log = open("../permutation.log_{}".format(thread_num),"w")
out_file = open("../promoter_go_terms/msy_term_enrichment_hypergeom_permutation_{}.table".format(thread_num),"w")
out_file.write("go_term\tp_val")

msy_enrichment_direction = []
msy_enrichment_p_vals = []
for perm in range(num_permut):
	start_time = time.time()
	rand_shuffle = range(len(gene_list))
	random.shuffle(rand_shuffle)
	reproduction = 0
	for j in range(len(go_terms)):
		#quesion 1 - what terms are enriched within the MSY region?
		msy_list = []
		msy_total = 0
		autosomal_list = []
		autosomal_total = 0
		for i in range(len(gene_list)):
			if i < 147:
				if abundance[i][j] > 0:
					msy_list.append(abundance[rand_shuffle[i]][j])
				#msy_total += tf_total[rand_shuffle[i]]
				msy_total += 1
			else:
				if abundance[i][j] > 0:
					autosomal_list.append(abundance[rand_shuffle[i]][j])
				#autosomal_total += tf_total[rand_shuffle[i]]
				autosomal_total += 1
		#len = presence, sum = abundance
		hpd = ss.hypergeom(autosomal_total,len(autosomal_list),msy_total)
		p_val = hpd.pmf(len(msy_list))
		out_file.write("\n{}\t{}".format(go_terms[j],p_val))
		if go_terms[j] == "GO:0000003":
			reproduction = p_val
	out_log.write("\nPermutation thread {}: replicate {} gave Reproduction = {}".format(thread_num,perm,reproduction))

out_log.close()
out_file.close()
		
