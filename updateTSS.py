import sys
import pandas as pd

xMap = dict()
with open('GeneLocus_Hg19.tsv', 'r') as fin:
	for line in fin:
		if line.startswith('#'):
			continue
		chunks = line.strip().split('\t')
		gene_id = chunks[0]
		if gene_id in xMap:
			print("ERROR: duplicate gene IDs.", gene_id, "\n")
			sys.exit(-1)
		else:
			xMap[gene_id] = line
res = []
with open('GREATv4.genes.hg19.tsv', 'r') as fin:
	for line in fin:
		if line.startswith('#'):
			continue
		chunks = line.strip().split('\t')
		CHROM = chunks[1]
		POS = chunks[2]
		STRAND = chunks[3]
		ids = (chunks[0]).split(',')
		if len(ids) == 1:
			tmp = ids[0]
			if tmp == 'NA':
				gene_id = chunks[4]
				gene_name = 'NA'
			else:
				gene_id = chunks[0]
				gene_name = chunks[4]
			if gene_id in xMap:
				block = (xMap[gene_id]).split('\t')
				ir = block[2]
				if gene_name == 'NA':
					gene_name = block[1]
			else:
				ir = ''
		else:
			length = []
			for tmp in ids:
				if tmp not in xMap:
					length.append(0)
					continue
				block = (xMap[tmp]).split('\t')[2]
				irTMP = block.replace('-',':').split(':')
				length.append(int(irTMP[2]) - int(irTMP[1]))
			flag = length.index(max(length))
			gene_id = ids[flag]
			block = (xMap[gene_id]).split('\t')
			gene_name = block[1]
			ir = block[2]
		text = [gene_id, CHROM, POS, STRAND, gene_name, ir]
		res.append(text)

df = pd.DataFrame(res, columns = ['gene_id', 'CHROM', 'POS', 'STRAND', 'gene_name', 'interval'])
df.to_csv('GeneTSS_Hg19.tsv', sep="\t", escapechar="\\", header=True, index=False, doublequote=False)
