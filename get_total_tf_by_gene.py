#!/usr/bin/python

import os

in_dir = os.listdir("../tf_list")
gene_list = []
tf_number = []
for in_file_name in in_dir:
        in_file = open("../tf_list/{}".format(in_file_name),"r")
	gene_list.append(in_file_name[:-4].strip(" "))
	line_count = -1
	for line in in_file:
		line_count += 1
	tf_number.append(line_count)
	in_file.close()

out_file = open("../promoter_go_terms/total_tf_by_gene.table","w")
out_file.write("gene_ID\ttf_count")
for i in range(len(gene_list)):
	out_file.write("\n{}\t{}".format(gene_list[i],tf_number[i]))
out_file.close()
