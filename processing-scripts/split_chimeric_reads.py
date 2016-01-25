#! /usr/bin/python

# Written by Simon Watson and Matt Cotten, Wellcome Trust Sanger Institute
# To use this script, change the following strings to point to your installs:
#
# 1) Location of QUASR Python modules folder (QUASR can be downloaded from https://github.com/simonjwatson/QUASR)
# 2) Location of BLAST database file stem
# 3) Location of BLASTN binary

import sys, os
sys.path.append('/Users/sw10/Dropbox/Sanger/QUASR/QUASR_v6.09/modules/') # 1)
import fastq
from Bio.Blast.Applications import NcbiblastnCommandline
from Bio.Blast import NCBIXML

def blast_reads(blast_string, reads, outfh):
	blast_db = '/Users/sw10/Dropbox/Sanger/blastdb/ebola/Zaire_ebolavirus_KM034562' # 2)
	blast_binary = '/Applications/ncbi-blast-2.2.29+/bin/blastn' # 3)
	xml_outfile = '/tmp/test.xml'
	evalue = 0.01 
	cline = NcbiblastnCommandline(cmd=blast_binary, out=xml_outfile, outfmt=5, query="-", db=blast_db, evalue=evalue, max_target_seqs=1, num_threads=1)
	stdout, stderr = cline(blast_string)

	with open(xml_outfile, 'r') as blast_handle:
		blast_records = NCBIXML.parse(blast_handle)
		for blast_record in blast_records:
			name = blast_record.query
			for alignment in blast_record.alignments:
				count = 1
				for hsp in alignment.hsps:
					seq = reads[name].sequence[hsp.query_start:hsp.query_end]
					qual = reads[name].quality[hsp.query_start:hsp.query_end]
					if hsp.sbjct_start > hsp.sbjct_end:
						tmp1 = [seq[i] for i in range(len(seq)-1,-1,-1)]
						seq = ''.join(tmp1)
						tmp2 = [qual[i] for i in range(len(qual)-1,-1,-1)]
						qual = ''.join(tmp2)
		
					outfh.write('@%s:%d\n%s\n+\n%s\n' % (name, count, seq, qual))
					count += 1
	os.remove(xml_outfile)


if len(sys.argv) != 4:
	print('split_chimeric_reads.py <infile.fq> <outfile.fq> <num_to_parse>')
	sys.exit(0)

infile = sys.argv[1]
outfile = sys.argv[2]
max_at_once = int(sys.argv[3])

counter = 0
reads = {}
blast_string = ''
with open(infile, 'r') as infh, open(outfile, 'w') as outfh:
	for header, sequence, quality in fastq.fastq_iterator(infh):
		blast_string += '>%s\n%s\n' % (header, sequence)
		reads[header] = fastq.FastqRecord(header, sequence, quality)
		counter += 1
		if counter == max_at_once:
			blast_reads(blast_string, reads, outfh)
			counter = 0
			reads = {}
			blast_string = ''
	blast_reads(blast_string, reads, outfh)
