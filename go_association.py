#!/usr/bin/python

from scipy import stats as ss
import numpy
import statsmodels.stats.multitest  as multitest
import time

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


in_file.close()

msy_enrichment_direction = []
msy_enrichment_p_vals = []
start_time = time.time()
for j in range(len(go_terms)):
	#quesion 1 - what terms are enriched within the MSY region?
	msy_list = []
	msy_total = 0
	autosomal_list = []
	autosomal_total = 0
	for i in range(len(gene_list)):
		if "|msy" in gene_list[i]:
			if abundance[i][j] > 0:
				msy_list.append(abundance[i][j])
			#msy_total += tf_total[i]
			msy_total += 1
		else:
			if abundance[i][j] > 0:
				autosomal_list.append(abundance[i][j])
			#autosomal_total += tf_total[i]
			autosomal_total += 1
	#len = presence, sum = abundance
	hpd = ss.hypergeom(autosomal_total,len(autosomal_list),msy_total)
	p_val = hpd.pmf(len(msy_list))
	msy_enrichment_p_vals.append(p_val)

print(time.time()-start_time)
out_file = open("../promoter_go_terms/msy_term_enrichment_hypergeom_noFWER.table","w")
out_file.write("go_term\tp_val")


for i in range(len(go_terms)):
	out_file.write("\n{}\t{}".format(go_terms[i],msy_enrichment_p_vals[i]))

out_file.close()
