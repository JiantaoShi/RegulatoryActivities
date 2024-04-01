import argparse
import time, sys
import pandas as pd
import numpy as np
import os, sys, pathlib

parser = argparse.ArgumentParser()
parser.add_argument('-A', "--annotation_file", type=str, default=None, required=True, help="Gene annotation file.")
parser.add_argument('-HWP', "--window_promoter", type=int, default=100, required=False, help="Half window size for promoters [100].")
parser.add_argument('-HWE', "--window_enhancer", type=int, default=10000, required=False, help="Half window size for enhancers [10000].")
parser.add_argument('-BN', "--bin_numbers", type=int, default=20, required=False, help="Number of bins for gene bodies [20].")
parser.add_argument('-C', "--cut_length", type=int, default=1000, required=False, help="Numober of bases to cut from gene bodies [1000].")
parser.add_argument('-L', "--min_length", type=int, default=1000, required=False, help="Minimal length of gene body [1000].")
parser.add_argument('-O', "--output_tag", type=str, default=None, required=True, help="Prefix for output files.")
args = parser.parse_args()

HWP = args.window_promoter
HWE = args.window_enhancer

# read gene annotation file
anno = []
with open(args.annotation_file, 'r') as fin:
	for line in fin:
		if line.startswith('gene_id'):
			continue
		anno.append(line.strip())

# make windows
def parseInterval(input):
	[CHR, START, END] = input.replace('-',':').split(':')
	return [int(START), int(END)]

res = []
for text in anno:
	[gene_id, CHROM, POS, STRAND, gene_name, gene_locus, reg_locus, promoter] = text.split('\t')
	POS = int(POS)
	gene_body = parseInterval(gene_locus)
	gene_regu = parseInterval(reg_locus)
	gene_prom = parseInterval(promoter)
	tag = gene_id + '_' + gene_name
	bed = [CHROM, POS - HWP, POS + HWP, tag + '_P_0']
	res.append(bed)
	if STRAND == '+':
		# promoter
		nu = round((POS - gene_prom[0])/HWP/2)
		for i in range(1, nu):
			center = POS - i*2*HWP
			bed = [CHROM, center - HWP, center + HWP, tag + '_PU_' + str(i)]
			res.append(bed)
		nd = round((gene_prom[1] - POS)/HWP/2)
		for i in range(1, nd):
			center = POS + i*2*HWP
			bed = [CHROM, center - HWP, center + HWP, tag + '_PD_' + str(i)]
			res.append(bed)
		# upstream distal
		nud = round((gene_prom[0] - gene_regu[0])/HWE/2)
		for i in range(0, nud):
			center = gene_prom[0] - i*2*HWE
			bed = [CHROM, center - HWE, center + HWE, tag + '_EU_' + str(i)]
			res.append(bed)
		# gene body
		gene_body[0] = gene_body[0] + args.cut_length
		bodySize = (gene_body[1] - gene_body[0])
		if bodySize > args.min_length:
			HWB = round(bodySize/args.bin_numbers/2)
			for i in range(1, args.bin_numbers + 1):
				center = gene_body[0] + i*2*HWB
				bed = [CHROM, center - HWB, center + HWB, tag + '_B_' + str(i)]
				res.append(bed)
		# downstream distal
		if gene_regu[1] - gene_body[1] > 2*HWE:
			ndd = round((gene_regu[1] - gene_body[1])/HWE/2)
			for i in range(1, ndd):
				center = gene_body[1] + 2*i*HWE
				bed = [CHROM, center - HWE, center + HWE, tag + '_ED_' + str(i)]
				res.append(bed)
	elif STRAND == '-':
		# promoter
		nd = round((POS - gene_prom[0])/HWP/2)
		for i in range(1, nd):
			center = POS - i*2*HWP
			bed = [CHROM, center - HWP, center + HWP, tag + '_PD_' + str(i)]
			res.append(bed)
		nu = round((gene_prom[1] - POS)/HWP/2)
		for i in range(1, nu):
			center = POS + 2*i*HWP
			bed = [CHROM, center - HWP, center + HWP, tag + '_PU_' + str(i)]
			res.append(bed)
		# upstream distal
		nud = round((gene_regu[1] - gene_prom[1])/HWE/2)
		for i in range(1, nud):
			center = gene_prom[1] + i*2*HWE
			bed = [CHROM, center - HWE, center + HWE, tag + '_EU_' + str(i)]
			res.append(bed)
		# gene body
		gene_body[1] = gene_body[1] - args.cut_length
		bodySize = (gene_body[1] - gene_body[0])
		if bodySize > args.min_length:
			HWB = round(bodySize/args.bin_numbers/2)
			for i in range(1, args.bin_numbers + 1):
				center = gene_body[1] - i*2*HWB
				bed = [CHROM, center - HWB, center + HWB, tag + '_B_' + str(i)]
				res.append(bed)
		# downstream distal
		if gene_body[0] - gene_regu[0] > 2*HWE:
			ndd = round((gene_body[0] - gene_regu[0])/HWE/2)
			for i in range(1, ndd):
				center = gene_body[0] - 2*i*HWE
				bed = [CHROM, center - HWE, center + HWE, tag + '_ED_' + str(i)]
				res.append(bed)

df = pd.DataFrame(res, columns = ['CHROM', 'start', 'end', 'name'])
df.to_csv(args.output_tag + '_windows.bed', sep="\t", escapechar="\\", header=False, index=False, doublequote=False)
