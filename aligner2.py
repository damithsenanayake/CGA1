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
    print "Usage: <reference file (fasta)> <read file 1 (fasta)> <read file 2 (fasta) > "
    sys.exit(0)

# Read the reference sequence and initiate the aligner
try:
    aligner = None

    for s in SeqIO.parse(sys.argv[1], "fasta"):
        aligner = Aligner1(s)
        #print aligner.refseq
        break  # Stop after the fist sequence in the reference
    out = []

    # read the first read file and find the mappings:

    R1s = []
    R2s = []

    for s in SeqIO.parse(sys.argv[2], "fasta"):
        R1s.append(s)
    for s in SeqIO.parse(sys.argv[3], "fasta"):
        R2s.append(s)

    for i in range(109):
        s1 = R1s[i]
        s2 = R2s[i]
        forward_out = ""
        reverse_out = ""
        # print (s.seq)
        ocount = 2
        alignment51 = pairwise2.align.localms(aligner.refseq, s1.seq, 2, -1, -1, -0.5)
        alignment52 = pairwise2.align.localms(aligner.refseq, s2.seq, 2, -1, -1, -0.5)
        alignment31 = pairwise2.align.localms(aligner.refseq.reverse_complement(), s1.seq, 2, -1, -1, -0.5)
        alignment32 = pairwise2.align.localms(aligner.refseq.reverse_complement(), s2.seq, 2, -1, -1, -0.5)
        score_1 = len(s1.seq._data) * 2
        score_2 = len(s2.seq._data) * 2
        count = 0
        r51 = 0
        r31 = 0
        r32 = 0
        r52 = 0
        c1 = False
        c2 = False
        strand = "+"
        position51 = len(aligner.refseq._data)
        position32 = len(aligner.refseq._data)
        ''' For the forward strand '''
        concordant = False

        concordant_alignments = []

        m31 = False
        m32 = False
        m51 = False
        m52 = False

        total_cons = 0;

        t51 = 1000000
        t52 = 1000000
        t31 = 1000000
        t32 = 1000000
        s1min = 100000
        s2min = 100000

        for a51 in alignment51:
            if a51[2] < score_1:
                continue
            s1min = min(s1min, a51[1]._data.index(s1.seq._data[0]))
        for a52 in alignment52:
            if a52[2] < score_2:
                continue
            s2min = min(s2min, a52[1]._data.index(s2.seq._data[0]))

        ordered = (s1min < s2min)

        #if ordered :
        for a51 in alignment51:
            if a51[2] == score_1:
                m51 = True
                t51 = a51[1]._data.index(s1.seq._data[0])

            for a32 in alignment32:
                if a32[2] == score_2:
                    m32 = True
                    t32 = a32[1]._data.index(s2.seq._data[0])

                    if (t51 + len(s1.seq._data)) < (len(aligner.refseq) - (t32 + len(s2.seq._data))):
                        concordant_alignments.append([t51+1, len(aligner.refseq) - (t32 + len(s2.seq._data))+1])

        #else :
        for a52 in alignment52:
            if a52[2] == score_2:
                m52 = True
                t52 = a52[1]._data.index(s2.seq._data[0])

            for a31 in alignment31:
                total_cons += 1
                if a31[2] == score_1:
                    m31 = True
                    t31 = a31[1]._data.index(s1.seq._data[0])

                    if (t52 + len(s2.seq._data)) < (len(aligner.refseq) - (t31 + len(s1.seq._data))):
                        concordant_alignments.append([t52+1, len(aligner.refseq) - (t31 + len(s1.seq._data))+1])

        #print concordant_alignments
        if len(concordant_alignments) > 0:
            count = len(concordant_alignments)
            position32 = len(aligner.refseq._data)
            position51 = len(aligner.refseq._data)
            for conc_alignment in concordant_alignments:
                if position32 > conc_alignment[1]:
                    position32 = conc_alignment[1]
                if position51 > conc_alignment[0]:
                    position51 = conc_alignment[0]

        if count > 0:
            out.append("r" + str(i + 1) + " " + aligner.refname + " " + str(position51) + " +" + " " + str(
                    count) + " " + aligner.refname + " " + str(position32) + " " + "-" + " " + str(count))
        else:
            out.append("r" + str(i + 1) + " * " + "0 " + "*" + " 0" + " * " + "0 " + "*" + " 0")
       # print total_cons
    with open('alignment2.txt', 'w') as f:
        for item in out:
            f.write("%s\n" % item)

except IOError as e:
    print "Could not read reference sequence file (see below)!"
    print e
    sys.exit(1)
