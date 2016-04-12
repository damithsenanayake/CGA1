import numpy as np

s = "alignment1.txt"
data = np.genfromtxt(s, dtype=[('read','S5'), ('ref','S5'), ('pos','i8'),('strand','S5'), ('count','i8')], delimiter=" ")

task1Text = "For task 1 in the assignment, the number of reads that exactly align with the reference is ;\n"

total_reads = np.unique(data['read']).shape[0]

task1Text += (str(total_reads) + "\n")

task1Text += "Out of them, the number of reads that align with the forward strand is: \n"

forward_count = len(np.where(data['strand']=='+')[0])
reverse_count = len(np.where(data['strand']=='-')[0])
task1Text += ("\t " + str(forward_count) + "\n")

task1Text += (" The number of reads that align with the reverse strand is :\n")
task1Text += ("\t" + str(reverse_count) + "\n")

task1Text += "\n The distribution of the alignments are as follows : \n" \
             "no_of_alignments\tread_count_with_no\n"

for i in np.unique(data['count']):
    task1Text += (str(i)+"\t\t\t\t\t"+str(len(np.where(data['count']==i)[0]))+":"+str(len(np.where(data['count']==i)[0])*100/total_reads)+"%\n")
print task1Text

f1 = open("analysis1.txt", "rw+")
f1.write(task1Text)

s = "alignment2.txt"
data = np.genfromtxt(s, dtype=[('read','S5'), ('ref','S5'), ('pos1','i8'),('strand1','S5'), ('count1','i8'), ('pos2','i8'),('strand2','S5'), ('count2','i8')], delimiter=" ")

task2Text = "For task 2 in the assignment, the number of reads that was considered is ;\n"

total_reads = np.unique(data['read']).shape[0]

task2Text += (str(total_reads) + "\n")

task2Text += "Out of them, the number of reads that have concordant alignments: \n"

forward_count = len(np.where(data['strand1']=='+')[0])
task2Text += ("\t " + str(forward_count) + "\n")


task2Text += "\n The distribution of number of concordant alignments is as follows : \n" \
             "no_of_alignments\tread_count_with_no\n"

for i in np.unique(data['count1']):
    task2Text += (str(i)+"\t\t\t\t\t"+str(len(np.where(data['count1']==i)[0]))+":"+str(len(np.where(data['count1']==i)[0])*100/total_reads)+"%\n")
f2 = open("analysis2.txt", "rw+")
f2.write(task2Text)


s = "alignment3.txt"
data = np.genfromtxt(s, dtype=[('read','S5'), ('ref','S5'), ('pos','i8'),('strand','S5'), ('count','i8'), ('cost', 'f8')], delimiter=" ")

task3Text = "For task 3 in the assignment, all the reads are expected to align with the reference. This is expected as we are penalizing\n" \
            "the inexact matchings in terms of a cost value. Therefore the total number of alignments in both forward and reverse strands is ;\n"

total_reads = np.unique(data['read']).shape[0]

task3Text += (str(total_reads) + "\n")


task3Text += "\n The distribution of the alignments are as follows : \n" \
             "no_of_alignments\tread_count_with_no\n"

for i in np.unique(data['count']):
    task3Text += (str(i)+"\t\t\t\t\t"+str(len(np.where(data['count']==i)[0]))+"\n")

task3Text += "\n The distribution of costs is as follows :\n" \
             "cost read_count_with_cost\n"

stats = []

for i in np.unique(data['cost']):
    task3Text += (str(i)+"\t\t\t\t\t"+str(len(np.where(data['cost']==i)[0]))+":"+str(len(np.where(data['count']==i)[0])*100/total_reads)+"%\n")
    for j in range (len(np.where(data['cost']==i)[0])):
        stats.append(i)

task3Text += "Additionally, the statistics(weighted by the number of alignments between a read and the reference) of the alignment costs are as below:\n"
task3Text += ("Highest observed cost = "+str(np.amax(data['cost']))+"\n")

task3Text += ("Median = " + str(np.median(np.array(stats)))+"\n")
task3Text += ("Mean = " + str(np.mean(np.array(stats))) + "\n")

task3Text += "\n It is also noteworthy that in some cases, a large number of alignments are reported with indel operations spanning a significant length\n" \
             "of the reference. In this case, we need to control the output space by introducing a gap penalty (for gap initiation and continuation)."

#calculate average


f3 = open("analysis3.txt", "rw+")
f3.write(task3Text)
