#!/usr/bin/python
from Bio.Seq import Seq
#load msy_contigs

in_msy_contigs = open("../../male_specific_contigs.bed","r")

msy_contigs = {}
autosomal_contigs = {}
for line in in_msy_contigs:
	msy_contigs[line.split("\t")[0]] = 1

in_msy_contigs.close()

#fill msy_contigs

in_genome = open("../../../reference_genome/male_genome/asm.contigs.fasta","r")

contig = ""
bases = ""
for line in in_genome:
	if line[0] == ">":
		if contig in msy_contigs:
			msy_contigs[contig] = bases
		else:
			autosomal_contigs[contig] = bases
		contig = line.strip(">").split(" ")[0]
		bases = ""
	else:
		bases += line.strip("\n")

in_genome.close()

#get genes in msy_contigs

in_gff = open("../../variants/combined/genes.gff","r")
out_file = open("../promoters/1kb_promoters.fa","w")
first = True
for line in in_gff:
	frags = line.split("\t")
	if frags[2] == "gene":
		if first:
			first = False
		else:
			out_file.write("\n")
		gene_name = frags[8].split(";")[0].split("=")[1]
		start = int(frags[3])
		end = int(frags[4])
		strand = frags[6]
		if frags[0] in msy_contigs:
			curr_seq = msy_contigs[frags[0]]
			gene_name += "|msy"
		else:
			curr_seq = autosomal_contigs[frags[0]]
			gene_name += "|autosomal"
		if strand == "+":
			if start < 1000:
				seq = curr_seq[0:start]
			else:
				seq = curr_seq[start-1000:start]
		else:
			if start+1000 > len(curr_seq):
				seq = Seq(curr_seq[start:])
			else:
				seq = Seq(curr_seq[start:start+1000])
			seq = seq.reverse_complement()
		out_file.write(">{}\n{}".format(gene_name,seq))

out_file.close()
in_gff.close()
