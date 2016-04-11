import numpy as np

class Wunscher(object):
    def __init__(self):
        self.scount = 1

    def match_string_build(self, A, B, seq1, seq2, grid, i, j):
        match = False
        if i > 0 and j > 0 and grid[i][j] == grid[i - 1][j - 1] + (0 if seq2[i - 1] == seq1[j - 1] else -1):
            A = seq2[i - 1] + A
            B = seq1[j - 1] + B
            i -= 1
            j -= 1
            match = True
        return A, B, i, j, match


    def delete_string_build(self, A, B, seq1, seq2, grid, i, j):
        match = False
        if j > 0 and grid[i][j] == grid[i][j - 1] - 1:
            B = seq1[j - 1] + B
            A = "_" + A
            j -= 1
            match = True
        return A, B, i, j, match


    def insert_string_build(self, A, B, seq1, seq2, grid, i, j):
        match = False
        if i > 0 and grid[i][j] == grid[i - 1][j] - 1:
            A = seq2[i - 1] + A
            B = "-" + B
            i -= 1
            match = True
        return A, B, i, j, match



    def backtrack_from_point(self, i, j, A, B, seq1, seq2, p):


        s = []
        visited ={}
        visited["0x0y"] =1
        s.append(str(i)+"x"+str(j)+"y")

        while len(s) > 0 :
            u = s.pop()
            if not u in visited.keys():
                visited[u]=1
                x = int(u.split("x")[0])
                y = int(u.split("x")[1].split("y")[0])

                diagParent = str(x-1)+"x"+str(y-1)+"y"

                if diagParent in p[u]:

                    s.append(diagParent)
                else:

                    for parent in p[u]:
                        #parent = str(parent[0])+"x"+str(parent[1])+"y"
                        if parent == "0x0y":# in visited.keys():
                            self.scount += 1
                            continue
                        s.append(parent)
            # else :
            #     self.scount +=1



        # if i == 0 and j == 0 :
        #     return self.scount + 1
        # parents = None
        # try:
        #     parents = p[str(i)+"x"+str(j)+"y"]
        # except KeyError:
        #     return 0
        # self.scount += (len(parents)-1)
        # for parent in parents:
        #
        #     self.backtrack_from_point(parent[0], parent[1], A, B, seq1, seq2,  p)
        #     #solutions.append([A,B])
        # return 0



    def BFS(self, graph,start,end):

        q = []
        validpaths =0
        temp_path = [start]

        q.append(temp_path)

        while not len(q) == 0:
            tmp_path = q[0]
            q = q[1:]
            last_node = tmp_path[len(tmp_path)-1]
            # print tmp_path
            if last_node == end:
                #print "VALID_PATH : ",tmp_path
                validpaths +=1
                #print validpaths

            try :
                for link_node in graph[last_node]:
                    if link_node not in tmp_path:
                        #new_path = []
                        new_path = tmp_path + [link_node]
                        q.append(new_path)
            except KeyError:
                continue

        return validpaths



    def needleman_align(self, seq1, seq2):
        seq1 = seq1[::-1]
        seq2 = seq2[::-1]

        parents = {}

        max_count_grid = np.zeros(shape=(len(seq2) + 1, len(seq1) + 1))

        grid = np.zeros(shape=(len(seq2) + 1, len(seq1) + 1))

        # initialize grid
        for i in range(grid.shape[0]):
            grid[i][0] = -1 * i

        for i in range(grid.shape[1]):
            grid[0][i] = -1 * i

            s1_out = "";
            s2_out = "";

        path = []
        for i in range(grid.shape[0]):
            # parents[str(i+1)+"x"+str(0)+"y"] =[[i, 0]]
            parents[str(i+1)+"x"+str(0)+"y"] =[str(i)+"x"+str(0)+"y"]

        for j in range(grid.shape[1]):
            # parents[str(0)+"x"+str(j+1)+"y"] = [[0, j]]
            parents[str(0)+"x"+str(j+1)+"y"] = [str(0)+"x"+str(j)+"y"]


        for i in range(1, grid.shape[0]):
            for j in range(1, grid.shape[1]):
                match = grid[i - 1][j - 1] + (0 if seq2[i - 1] == seq1[j - 1] else -1)
                delete = grid[i - 1][j] - 1
                insert = grid[i][j - 1] - 1

                grid[i][j] = max(match, delete, insert)
                parents[str(i)+"x"+str(j)+"y"] = []
                maxcount = 0
                if match == grid[i][j] :
                    maxcount += 1
                    # parents[str(i)+"x"+str(j)+"y"].append([i-1, j-1])
                    parents[str(i)+"x"+str(j)+"y"].append(str(i-1)+"x"+str(j-1)+"y")

                if delete == grid[i][j] :
                    maxcount += 1
                    # parents[str(i)+"x"+str(j)+"y"].append([i-1, j])
                    parents[str(i)+"x"+str(j)+"y"].append(str(i-1)+"x"+str(j)+"y")

                if insert == grid[i][j] :
                    maxcount += 1
                    # parents[str(i)+"x"+str(j)+"y"].append([i, j-1])
                    parents[str(i)+"x"+str(j)+"y"].append(str(i)+"x"+str(j-1)+"y")

                max_count_grid[i][j] = maxcount


        A = ""
        B = ""
        i = len(seq2)
        j = len(seq1)

        A, B = self.backtrack_for_optimal(A, B, grid, i, j, seq1, seq2)
        #print parents
        #print grid
        self.backtrack_from_point(i, j, "", "", seq2, seq1, parents)
        #self.scount = self.BFS(parents,str(len(seq2))+"x"+str(len(seq1))+"y", "0x0y")
        #self.count_paths(parents,str(len(seq2))+"x"+str(len(seq1))+"y", "0x0y")
        #print B[::-1]
        #print A[::-1]
        return B[::-1], A[::-1],self.scount, grid[len(seq2)][len(seq1)]



        #print self.scount
    '''public AllPaths(Graph G, String s, String t) {
        enumerate(G, s, t);
    }

    // use DFS
    private void enumerate(Graph G, String v, String t) {

        // add node v to current path from s
        path.push(v);
        onPath.add(v);

        // found path from s to t - currently prints in reverse order because of stack
        if (v.equals(t))
            StdOut.println(path);

        // consider all neighbors that would continue path with repeating a node
        else {
            for (String w : G.adjacentTo(v)) {
                if (!onPath.contains(w)) enumerate(G, w, t);
            }
        }

        // done exploring from v, so remove from path
        path.pop();
        onPath.delete(v);
    }'''

    def count_paths(self, graph, source, destination):

        path = []
        onpath = []

        path.append(source)
        onpath.append(source)

        if source == destination:
            self.scount += 1
        else :
            for w in graph[source]:
                if not w in onpath :
                    self.count_paths(graph, w, destination)
        path.pop()
        onpath.remove(source)

    def backtrack_for_optimal(self, A, B, grid, i, j, seq1, seq2):
        while i > 0 or j > 0:
            # self.scount += (mc_grid[i][j] - 1 )
            A, B, i, j, match = self.match_string_build(A, B, seq1, seq2, grid, i, j)

            delete = False;
            if not match:
                A, B, i, j, delete = self.delete_string_build(A, B, seq1, seq2, grid, i, j)

            if not delete:
                A, B, i, j, insert = self.insert_string_build(A, B, seq1, seq2, grid, i, j)

        return A, B


#Wunscher().needleman_align("TCGTAGATTTTTGACCCCCGA", "GA")

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
        print aligner.refseq
        break  # Stop after the fist sequence in the reference
    out = []
    for s in SeqIO.parse(sys.argv[2], "fasta"):
        forward_out = ""
        reverse_out = ""
        print (s.seq)
        ocount=2
        A, B, count, score = Wunscher().needleman_align(aligner.refseq._data, s.seq._data)
        strand = "+"
        position = B.index(s.seq._data[0])
        ''' For the forward strand '''
        sim = len(aligner.refseq) + score
        cost = len(s) - sim

        forward_out = s.id +  " " + aligner.refname + " "+ str(position+1) + " " + strand + " " +str(count) + " " + str(cost)
        ''' For the Reverse Strand '''
        position = len(aligner.refseq._data)
        count = 0
        strand = "-"
        A, B, count, score = Wunscher().needleman_align(aligner.refseq.reverse_complement().complement()._data, s.seq._data)

        sim = len(aligner.refseq) + score
        cost = len(s) - sim

        position = B.index(s.seq.reverse_complement().complement()._data[0])

        reverse_out = s.id +  " " + aligner.refname + " "+ str(position+1) + " " + strand + " " +str(count) + " " + str(cost)

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

