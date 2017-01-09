


import os, sys

def read_meth_file(f, n):
	"""
		Given the methylation file and the chromosome number(str).
		Return a dictionary: 
		key is methylation position(int);
		value is a list of methylation type, num_C(int), depth(int).
	"""
	res = {}
	for line in f:
		line = line.split()
		if line[0] == n:
			pos = int(line[1]) # pos is an integer
			res[pos] = [line[3],int(line[4]),int(line[5])] # e.g.['CHG',3,4]
	return res


def read_gff_file(f, n):
	"""
		Given the GFF file and the chromosome number(str).
		Return a dictionary:
		key is the gene name;
		value is tuple, '+': first<second,'-': first>second.
	"""
	res = {}
	for line in f:
		line = line.split()
		if line[0][:2] == n:
			g_name = line[2]
			sp,ep = int(line[3]),int(line[4])
			strand = line[6]
			if strand == '+':
				res[g_name] = (sp, ep)
			else:
				res[g_name] = (ep, sp) # if on '-' strand, starting pos is bigger
	return res

def construct_bins(pos,bs,bn):
	"""
		Given the postion of the gene and the strand.
		Return the starting postion of each bin.
	"""
	binfor5, binfor3 = [],[]
	if pos[0] < pos[1]:
		for b in range(bn, -bn, -1):
			binfor5.append(pos[0] - bs*b)
			binfor3.append(pos[1] + 1 - bs*b)
	else:
		for b in range(bn, -bn, -1):
			binfor5.append(pos[0] + bs*b)
			binfor3.append(pos[1] - 1 + bs*b)

	return binfor5,binfor3

def write_meth_each_bin(meths, bins, fl, bs):
	"""
		Given the list of bin positions.
		Write the files in filehandle list.
		No return value.
	"""
	if bins[0] > bins[-1]:
		label = '-'
	else:
		label = '+'
	for i in range(len(bins)):
		if label == '+':
			each_bin = [bins[i]+j for j in range(bs)]
		else:
			each_bin = [bins[i]-j for j in range(bs)]
		total = [0,0,0,0] # [totalMeth,totalCG,totalCHG,totalCHH]
		cnum = [0.,0.,0.,0.] # [Cnum,CGnum,CHGnum,CHHnum]
		for eb in each_bin:
			if eb not in meths: # not having methylation or beyond chromosome
				continue
			total[0] += meths[eb][2]
			cnum[0] += meths[eb][1]

			if meths[eb][0] == 'CG':
				cnum[1] += meths[eb][1]
				total[1] += meths[eb][2]
			elif meths[eb][0] == 'CHG':
				cnum[2] += meths[eb][1]
				total[2] += meths[eb][2]
			else:
				cnum[3] += meths[eb][1]
				total[3] += meths[eb][2]
		
		for x in [0,1,2,3]:
			f = fl[x]
			if total[x] == 0:
				f.write('\t%d' % total[x])
			else:
				f.write('\t%2.4f' % (100.*cnum[x]/total[x]))
	for f in fl:
		f.write('\n')
	return None




def write_meth_level_for_every_gene(genes,meths,fl,bs,bn):
	for gene in genes:
		for f in fl:
			f.write('%s' % gene)
		locs = genes[gene]
		
		binfor5, binfor3 = construct_bins(locs,bs,bn)
		# print gene,locs,binfor5,binfor3
		write_meth_each_bin(meths,binfor5,fl[:4],bs)
		write_meth_each_bin(meths,binfor3,fl[4:],bs)
	return None


if len(sys.argv) != 7:
	raise SyntaxError('Input should include methylation file, gff file, chromosome number, bin size, bin number and the mutant out pre!')

input_name,gff_name,chromosome_n,bin_size,bin_number,output_pre = sys.argv[1],sys.argv[2],sys.argv[3],sys.argv[4],sys.argv[5],sys.argv[6]
bin_size = int(bin_size)
bin_number = int(bin_number)
chromosome_n = int(chromosome_n)


file_h = [] # list containing the output filehandles
for i in ['5','3']:
	for j in ['c','cg','chg','chh']:
		name = output_pre+'_'+j+'_'+i+'.txt' # e.g. 'mutant_c_3.txt','mutant_chh_5.txt'
		f_h = open(name, 'w')
		file_h.append(f_h)


for i in range(chromosome_n+1): # include chr00 as the unknown
	input_f = open(input_name)
	input_f.readline() # discard header
	gff_f = open(gff_name)
	chr_number = '0'+str(i)
	chr_number = chr_number[-2:] # 00,01,02...11,12
	meth_dic = read_meth_file(input_f, chr_number)
	gene_dic = read_gff_file(gff_f, chr_number)
	write_meth_level_for_every_gene(gene_dic,meth_dic,file_h,bin_size,bin_number)
	input_f.close()
	gff_f.close()


for f in file_h:
	f.close()









