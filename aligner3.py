import numpy as np


class SWAligner(object):
    def __init__(self):
        self.scount = 0



    def modified_sw_align(self, seq1, seq2):
        seq1 = seq1[::-1]
        seq2 = seq2[::-1]
        grid = np.zeros(shape=(len(seq2) + 1, len(seq1) + 1))


        for i in range(1, grid.shape[0]):
            for j in range(1, grid.shape[1]):
                match = grid[i - 1][j - 1] + (1 if seq2[i - 1] == seq1[j - 1] else 0)
                delete = grid[i - 1][j] - 1
                insert = grid[i][j - 1] - 1

                grid[i][j] = max(match, delete, insert)


        position, count = self.backtrack_and_count(grid)

        return count,position,np.amax(grid)



    def backtrack_and_count(self, grid):

        j = grid.shape[0]-1
        first_position = 0;
        maximum = np.amax(grid)
        count = 0
        first_found = False

        for i in xrange(grid.shape[1]-1, 0, -1):
            if grid[j][i] >= grid[j-1][i-1] and grid[j][i] == maximum :
                count += 1
                if not first_found:
                    first_position = grid.shape[1]-i
                    first_found = True


        return first_position, count








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
        # print alignment
        self.alignments.append(alignment)
        return self.alignments


# Check the command line arguments
sys.setrecursionlimit(10000)
if len(sys.argv) < 3:
    print "Usage: <reference file (fasta)> <read file (fasta)> "
    sys.exit(0)

# Read the reference sequence and initiate the aligner
try:
    aligner = None

    for s in SeqIO.parse(sys.argv[1], "fasta"):
        aligner = Aligner1(s)
        #print aligner.refseq
        break  # Stop after the fist sequence in the reference
    out = []
    for s in SeqIO.parse(sys.argv[2], "fasta"):
        forward_out = ""
        reverse_out = ""
        #print (s.seq)
        ocount = 2
        swaligner = SWAligner()
        count,position, score = swaligner.modified_sw_align(aligner.refseq._data, s.seq._data)
        strand = "+"
        ''' For the forward strand '''

        forward_out = s.id + " " + aligner.refname + " " + str(position) + " " + strand + " " + str(
                count) + " " + str(len(s)-score)
        ''' For the Reverse Strand '''

        swaligner = SWAligner()
        position = len(aligner.refseq._data)
        count = 0
        strand = "-"
        count,position, score = swaligner.modified_sw_align(aligner.refseq.reverse_complement().complement()._data,
                                                                  s.seq._data)


        reverse_out = s.id + " " + aligner.refname + " " + str(position) + " " + strand + " " + str(
                count) + " " + str(len(s)-score)

        if ocount != 0:
            if forward_out != None:
                out.append(forward_out)
            if reverse_out != None:
                out.append(reverse_out)
        else:
            out.append(s.id + "* " + "0 " + "*" + " 0")
    with open('alignment3.txt', 'w') as f:
        for item in out:
            f.write("%s\n" % item)
except IOError as e:
    print "Could not read reference sequence file (see below)!"
    print e
    sys.exit(1)
