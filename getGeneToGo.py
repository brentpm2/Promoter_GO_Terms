#!/usr/bin/python

in_file = open("../promoter_annotation/go-to-gene.map","r")

mapper = {}

for line in in_file:
	line = line.strip("\n")
	frags = line.split(" ")
	is_name = True
	name = ""
	for frag in frags:
		if is_name:
			is_name = False
			name = frag
		else:
			if frag in mapper:
				mapper[frag].append(name)
			else:
				mapper[frag] = [name]
in_file.close()

out_file = open("../promoter_annotation/gene-to-go.map","w")

first_line = True
for key in mapper:
	if first_line:
		first_line = False
	else:
		out_file.write("\n")
	out_file.write(key)
	for element in mapper[key]:
		out_file.write(" ")
		out_file.write(element)
out_file.close()
