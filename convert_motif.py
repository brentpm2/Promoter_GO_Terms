#!/usr/bin/python

import os

in_mapping = open("../promoter_annotation/ID_mapping_Arabidopsis_thaliana.txt","r")

converter = {}

first = True
for line in in_mapping:
	if first:
		first = False
	else:
		frags = line.split("\t")
		if frags[1] != "NA":
			converter[frags[0]] = frags[1]
in_mapping.close()

#in_file = open("../promoter_motifs/ maleWH.00g000010|autosomal .gff","r")
#out_file = open("../tf_list/maleWH.00g000010|autosomal.txt","w")

in_dir = os.listdir("../promoter_motifs/")

first_in = True
first_out = True
for in_file_name in in_dir:
        in_file = open("../promoter_motifs/{}".format(in_file_name),"r")
        out_file = open("../tf_list/{}".format(in_file_name),"w")
        for line in in_file:
	        if first_in:
		        first_in = False
	        else:
		        frags = line.split("\t")
		        if frags[0] in converter:
			        if first_out:
				        first_out=False
			        else:
				        out_file.write("\n")
			        out_file.write(converter[frags[0]])

        in_file.close()
        out_file.close()
