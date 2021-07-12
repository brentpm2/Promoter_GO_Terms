#!/usr/bin/python


in_file = open("../promoter_go_terms/GO_to_description.map","r")
desc_map = {}
for line in in_file:
	line = line.strip("\n")
	frags = line.split("\t")
	desc_map[frags[0]] = frags[1]
in_file.close()

in_file = open("../promoter_go_terms/msy_term_enrichment_hypergeom_FWER_permutation.sig.table","r")
out_file = open("../promoter_go_terms/msy_term_enrichment_hypergeom_FWER_permutation.presence.description.sig.table","w")
out_file.write("GO_term\tadj.p_val\tdescription")
first = True
for line in in_file:
	if first:
		first = False
	else:
		line = line.strip("\n")
		frags = line.split("\t")
		out_file.write("\n{}\t{}".format(line,desc_map[frags[0]]))
in_file.close()
out_file.close()
