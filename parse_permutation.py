#!/usr/bin/python

import os
import statsmodels.stats.multitest  as multitest

in_dir = os.listdir("../promoter_go_terms")

go_terms = {}
first = True

for in_file_name in in_dir:
	if "msy_term_enrichment_hypergeom_permutation_" in in_file_name:
		in_file = open("../promoter_go_terms/{}".format(in_file_name),"r")
		first = True
		for line in in_file:
			if first:
				first = False
			else:
				line = line.strip("\n")
				frags = line.split("\t")
				if frags[1] != "nan":
					if frags[0] in go_terms:
						go_terms[frags[0]].append(float(frags[1]))
					else:
						go_terms[frags[0]] = [float(frags[1])]
		in_file.close()



in_file = open("../promoter_go_terms/msy_term_enrichment_hypergeom_noFWER.table","r")

present_terms = []
present_terms_out = open("../promoter_go_terms/present_go_terms.list","w")

first = True
first_out = True
p_vals = []
go_term_ID = []
for line in in_file:
	if first:
		first = False
	else:
		line = line.strip("\n")
		frags = line.split("\t")
		if frags[0] in go_terms:
			p_vals.append(sum(i <= float(frags[1]) for i in go_terms[frags[0]]) * 1.0 / len(go_terms[frags[0]]))
			go_term_ID.append(frags[0])
			if first_out:
				first_out = False
			else:
				present_terms_out.write("\n")
			present_terms_out.write(frags[0])
			
		else:
			print("Error, missing term: {}".format(frags[0]))
in_file.close()
print(len(p_vals),p_vals[0])
fdr_correction = multitest.fdrcorrection(p_vals)
out_file = open("../promoter_go_terms/msy_term_enrichment_hypergeom_FWER_permutation.table","w")
out_file.write("go_term\tp_val")
out_file_sig = open("../promoter_go_terms/msy_term_enrichment_hypergeom_FWER_permutation.sig.table","w")
out_file_sig.write("go_term\tp_val")
num_sig = 0

for i in range(len(fdr_correction[1])):
	out_file.write("\n{}\t{}".format(go_term_ID[i],fdr_correction[1][i]))
	if fdr_correction[1][i] < 0.05:
		out_file_sig.write("\n{}\t{}".format(go_term_ID[i],fdr_correction[1][i]))
		num_sig += 1

present_terms_out.close()
out_file.close()
out_file_sig.close()
print(num_sig)
