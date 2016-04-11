#!/usr/bin/env python
import sys
from Bio import SeqIO, pairwise2


# Support Class to create alignments and check relationships between alignments
from Bio.pairwise2 import format_alignment


class Alignment:
    def __init__(self, chr, pos, strand):
        self.chr = chr
        self.pos = pos
        self.strand = strand

    def is_concordant_with(self, other):
        return False

    def __repr__(self):
        return self.chr + "\t" + str(self.pos) + "\t" + self.strand


# Support class to execute the handiwork of read alignment to a reference sequence
class Aligner1:
    def __init__(self, ref):
        self.refname = ref.id
        self.refseq = ref.seq
        self.rrefseq = ref.seq.reverse_complement()
        self.alignments = []

    def align(self, read):
        alignment = pairwise2.align.globalxx(self.refseq, read)
        print alignment
        self.alignments.append(alignment)
        return self.alignments


# Check the command line arguments
if len(sys.argv) < 3:
    print "Usage: <reference file (fasta)> <read file (fasta)> "
    sys.exit(0)

# Read the reference sequence and initiate the aligner
try:
    aligner = None

    for s in SeqIO.parse(sys.argv[1], "fasta"):
        aligner = Aligner1(s)
        print aligner.refseq
        break  # Stop after the fist sequence in the reference
    out = []
    for s in SeqIO.parse(sys.argv[2], "fasta"):
        forward_out = ""
        reverse_out = ""
        #print (s.seq)
        ocount=2
        alignment = pairwise2.align.localms(aligner.refseq, s.seq, 2, -1, -1, -0.5)
        count = 0
        strand = "+"
        position = len(aligner.refseq._data)
        ''' For the forward strand '''
        for a in alignment:
            if a[2] < 12:
                continue
            t = a[1]._data.index(s.seq._data[0])
            if position > t :
                position = t
            count += 1
            #print(format_alignment(*a))
        if count == 0:
            ocount-=1
            position = 0
            forward_out = None
        else :
            forward_out = s.id +  " " + aligner.refname + " "+ str(position+1) + " " + strand + " " +str(count)
        ''' For the Reverse Strand '''
        position = len(aligner.refseq._data)
        count = 0
        strand = "-"
        alignment = pairwise2.align.localms( aligner.refseq, s.seq.reverse_complement().complement(), 2, -1, -1, -0.5 )
        for a in alignment :
            if a[2] < 12:
                continue
            t = a[1]._data.index(s.seq.reverse_complement().complement()._data[0])
            if position > t :
                position = t
            count += 1
            #print(format_alignment(*a))
        if count == 0:
            ocount-=1
            position = 0
            reverse_out = None
        else :
            reverse_out = s.id + " " + aligner.refname + " "+ str(position + 1) + " " + strand + " " +str(count)

        if ocount != 0 :
            if forward_out != None :
                out.append(forward_out)
            if reverse_out != None :
                out.append(reverse_out)
        else :
            out.append(s.id + "* " + "0 " + "*" + " 0")
    with open('alignment1.txt','w') as f :
        for item in out:
            f.write("%s\n" % item)
except IOError as e:
    print "Could not read reference sequence file (see below)!"
    print e
    sys.exit(1)
