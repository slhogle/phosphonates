#!/usr/bin/env python
import sys
import random
import pandas as pd
import numpy as np
from Bio import SeqIO


def id_reader(family):
	with open('/nobackup1/shogle/seq_databases/diamond/MARMICRODB-dedup/seq_families/pro-syn-pel-kaks10/'+family+'/'+family+'.ids') as f_in:
		return [line.rstrip() for line in f_in]

def faa_reader(family):
	myfa = '/nobackup1/shogle/seq_databases/diamond/MARMICRODB-dedup/seq_families/pro-syn-pel-kaks10/'+family+'/'+family+'.faa'
	lens = [len(record.seq) for record in SeqIO.parse(myfa, "fasta")]
	return np.median(lens)*3

def checkEqual(lst):
	return not lst or lst.count(lst[0]) == len(lst)

def parse_diamond(m8):
	""" Yield formatted record from diamond custom field m8 file """
	formats = {0:str, 1:str, 2:float, 3:float, 4:float, 5:float, 6:float, 7:float, 8:float, 9:float}
	fields = {0:'read',1:'dbtarget',2:'bitscore',3:'evalue',4:'pident',5:'qcovhsp',6:'dbtargetlength_aa',7:'readlength_nt',8:'alnlength_aa'}
	with open(m8) as f_in:
		for line in f_in:
			yield dict([(fields[index], formats[index](value)) for index, value in enumerate(line.rstrip().split())])

def find_top_hits(diamond_file):
	x = {}
	for y in parse_diamond(diamond_file):
		if y['qcovhsp'] < 50: # does not meet marker percentage alignment cutoff
			continue
		elif y['read'] not in x: # record aln
			x[y['read']] = [y]
		elif x[y['read']][0]['bitscore'] == y['bitscore']: # add aln if score is the best score of a group
			x[y['read']] += [y]
		elif x[y['read']][0]['bitscore'] < y['bitscore']: # update aln if new one is better
			x[y['read']] = [y]
	return x

def threshold4hits(all_best_hits, markers, threshold):
	"""Filter best hit collections for each read so that at least X threshold of top scoring alignments hit a marker"""
	"""Return only besthit alignments matching a marker family"""
	w = {}
	for x,y in all_best_hits.iteritems():
		q = [z for z in y if z['dbtarget'] in markers.keys()]
		if len(q)/len(y) > threshold:
			w[x] = q
	return w

def unambiguous_hits(all_best_hits, markers):
	ua = dict([(_,[]) for _ in markers.values()])
	for x,y in all_best_hits.iteritems():
		f = [markers[z['dbtarget']] for z in y]
		if checkEqual(f):
			y[0]['dbtargetlength_aa'] = np.median([z['dbtargetlength_aa'] for z in y])
			ua[f[0]].append(y[0])
	return ua

def ambiguous_hits(all_best_hits, unambig_hits, markers):
	""" Probabalistically assign ambiguously mapped reads """
	final_hits = unambig_hits.copy()
	for x,y in all_best_hits.iteritems():
		f = [markers[z['dbtarget']] for z in y]
		if not checkEqual(f): # more than one hit
			if checkEqual([i[4:] for i in list(set(f))]): # but the protein family of the hits is the same
				q = list(set(f))
				counts = [len(final_hits[p]) for p in q]
				if sum(counts) == 0:
					assignedmarker = random.sample(q, 1)[0]
				else:
					probs = [float(count)/sum(counts) for count in counts]
					assignedmarker = np.random.choice(q, 1, p=probs)[0]
				y[0]['dbtargetlength_aa'] = np.median([z['dbtargetlength_aa'] for z in y])
				final_hits[assignedmarker].append(y[0])
	return final_hits

def report_rpk_exact(alns):
	t=[]
	for x in alns:
		t.append(1/(x['dbtargetlength_aa']*3/1000))
	return sum(t)

def report_rpk_median(alns, fam_l):
	t=[]
	for x in alns:
		t.append( 1/(fam_l/1000))
	return sum(t)

######################################################
######################## MAIN ########################
######################################################
station=sys.argv[1]
meta=sys.argv[2]

base='/nobackup1b/users/shogle/projects/pro-syn-sar-core-counts/raw_output/'+meta
reports='/nobackup1/shogle/projects/pro-syn-sar-core-counts/reports/'+meta

#stations=["S{:04d}".format(x) for x in range(1, 630)]
families=['pel-PX000099', 'pel-PX000249', 'pel-PX000355', 'pel-PX000369', 'pel-PX000374', 'pel-PX000396', 'pel-PX000397', 'pel-PX000612', 'pel-PX000669', 'pel-PX000683', 'pro-PX000049', 'pro-PX000092', 'pro-PX000116', 'pro-PX000176', 'pro-PX000234', 'pro-PX000373', 'pro-PX000453', 'pro-PX000495', 'pro-PX000732', 'pro-PX000753', 'syn-PX000049', 'syn-PX000092', 'syn-PX000116', 'syn-PX000176', 'syn-PX000234', 'syn-PX000373', 'syn-PX000453', 'syn-PX000495', 'syn-PX000732', 'syn-PX000753']

fam_d = {}
for f in families:
	for i in id_reader(f):
		fam_d[i]=f

fam_l = {}
for f in families:
	fam_l[f] = faa_reader(f)

#report = dict([(_,{}) for _ in stations])
report = {station: {}}

all_best_hits = find_top_hits(base+'/'+str(station)+'-ALL-RECIPR.m6')
filtered_best_hits = threshold4hits(all_best_hits, fam_d, 0.75)
uahits = unambiguous_hits(filtered_best_hits, fam_d)
final = ambiguous_hits(filtered_best_hits, uahits, fam_d)
for f,h in final.iteritems():
	report[station].update({f: report_rpk_median(h, fam_l[f])})

(pd.DataFrame.from_dict(report, orient="index").fillna(value=0))[families].to_csv(reports+"/"+station+".tsv", header=False, sep='\t', index=True)

