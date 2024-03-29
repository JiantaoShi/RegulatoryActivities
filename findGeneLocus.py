import gzip
import csv
import pandas as pd

kchr = list(map(str, range(1,23)))
kchr.extend(['MT', 'X', 'Y'])

res = []
with gzip.open('Homo_sapiens.GRCh37.87.gtf.gz','rt') as fin:        
    for line in fin:        
        if line.startswith('#'):
        	continue
        chunks = line.strip().split('\t')
        if chunks[2] != 'gene':
        	continue
        strand = chunks[6]
        if chunks[0] not in kchr:
        	continue
        chr = 'chr' + chunks[0]
        if chr == 'chrMT':
        	chr = 'chrM'
        ir = chr + ':' + chunks[3] + '-' + chunks[4]
        anno = (chunks[8]).split(';')
        gene_id = None
        gene_name = None
        for block in anno:
        	md = block.replace('"', '').lstrip().split(' ')
        	if md[0] == 'gene_id':
        		gene_id = md[1]
        	elif md[0] == 'gene_name':
        		gene_name = md[1]
        if gene_id is not None and gene_name is not None:
        	text = [gene_id, gene_name, ir, strand]
        	res.append(text)

df = pd.DataFrame(res, columns = ['gene_id', 'gene_name', 'interval', 'strand'])
df.to_csv('GeneLocus_Hg19.tsv', sep="\t", escapechar="\\", header=True, index=False, doublequote=False)
