sequence = ('.txt')
#sequence = open('seq.txt',mode='rb')
from DNA_Reverse_Compliment import revCompIterative
Found_ORFs = open('./Predicted_ORFs.txt',mode='wb')
Longer_ORFs = open('./Longer_Predicted_ORFs.txt',mode='wb')
Larger_ORFs = open('./Larger_Predicted_ORFs.txt',mode='wb')
import re
import collections
import sys

def find_all(sequence, subsequence):
    ''' Returns a list of indexes within sequence that are the start of subsequence'''
    start = 0
    idxs = []
    next_idx = sequence.find(subsequence, start)

    while next_idx != -1:
        idxs.append(next_idx)
        start = next_idx + 1  # Move past this on the next time around
        next_idx = sequence.find(subsequence, start)

    return idxs




#fname = file(sys.argv[1])  # Read in from the first command-line argument

Contigs = []

#seq = sequence.read()
Contig = ""

#with open(inSequence, 'rb') as line:
# with open(sequence) as f:
#     lines = f.readlines()
#
#     for line in lines:
#         if line is lines[-1]:
#             Contigs.append(Contig)
#             Contig = ""
#
#         line = line.replace("\n", "")
#         if ">" not in line:
#             Contig += str(line)
#         elif ">" in line:
#             Contigs.append(Contig)
#             Contig = ""
#         elif line is None:
#             print "end"
#             Contigs.append(Contig)
#             Contig = ""
#
# Genome = str(Contigs)
# Genome = filter(str.isalnum, Genome)
Genome = ""
with open (sequence, 'rb') as Gfile:
    for line in Gfile:
        line = line.replace("\n","")
        if ">" not in line:
            Genome += str(line)

print Genome
Genome_rev = revCompIterative(Genome)

print "stops"
print Genome_rev

stop_amber = find_all(Genome, 'TAG')
stop_ochre = find_all(Genome, 'TAA')
stop_umber = find_all(Genome, 'TGA')
stops = stop_amber + stop_ochre + stop_umber
stops.sort()


Frame_1 = collections.OrderedDict()
Frame_2 = collections.OrderedDict()
Frame_3 = collections.OrderedDict()

Seen_Stops = []
counter = 0
Stop_Stop = collections.OrderedDict()

Lengths = []

for stop in stops:
    Found = False
    #if stop not in Seen_Stops:
    Seen_Stops.append(stop)
    if stop != stops[-1]:
        for next_stop in stops[counter+1:]:
            diff = abs(stop-next_stop)
            Lengths.append(diff)

            if diff % 3 == 0:
                if next_stop not in Seen_Stops:

                    #if diff > 500:
                    frame = stop % 3
                    if frame == 0:
                        Frame_1.update({stop:next_stop})
                    if frame == 1:
                        Frame_2.update({stop:next_stop})
                    if frame == 2:
                        Frame_3.update({stop:next_stop})

                    Found = True
                    Seen_Stops.append(next_stop)
                    break
    # if Found == Fals:
    #     Stop_Stop

    counter +=1

print max(Lengths)
#################################
Found_ORFs.write("Frame 1 \n")
for k, v in Frame_1.iteritems():
    Found_ORFs.write(str(k) + "\t" + "+" + "\t" + str(v) + "\n")
Found_ORFs.write("Frame 2 \n")
for k, v in Frame_2.iteritems():
    Found_ORFs.write(str(k) + "\t" + "+" + "\t" + str(v) + "\n")
Found_ORFs.write("Frame 3 \n")
for k, v in Frame_3.iteritems():
    Found_ORFs.write(str(k) + "\t" + "+" + "\t" + str(v) + "\n")
#################################
Longer_ORFs.write("Frame 1 \n")
for k, v in Frame_1.iteritems():
    if abs(k-v) > 200:
        Longer_ORFs.write(str(k) + "\t" + "+" + "\t" + str(v) + "\n")
Longer_ORFs.write("Frame 2 \n")
for k, v in Frame_2.iteritems():
    if abs(k - v) > 200:
        Longer_ORFs.write(str(k) + "\t" + "+" + "\t" + str(v) + "\n")
Longer_ORFs.write("Frame 3 \n")
for k, v in Frame_3.iteritems():
    if abs(k-v) > 200:
        Longer_ORFs.write(str(k) + "\t" + "+" + "\t" + str(v) + "\n")
##################################
Larger_ORFs.write("Frame 1 \n")
for k, v in Frame_1.iteritems():
    if abs(k - v) >= 350:
        Larger_ORFs.write(str(k) + "\t" + "+" + "\t" + str(v) + "\n")
Larger_ORFs.write("Frame 2 \n")
for k, v in Frame_2.iteritems():
    if abs(k - v) >= 350:
        Larger_ORFs.write(str(k) + "\t" + "+" + "\t" + str(v) + "\n")
Larger_ORFs.write("Frame 3 \n")
for k, v in Frame_3.iteritems():
    if abs(k - v) >= 350:
        Larger_ORFs.write(str(k) + "\t" + "+" + "\t" + str(v) + "\n")

###### Reversed
Genome_rev = revCompIterative(Genome)

print Genome_rev

stop_amber = find_all(Genome_rev, 'TAG')
stop_ochre = find_all(Genome_rev, 'TAA')
stop_umber = find_all(Genome_rev, 'TGA')
stops = stop_amber + stop_ochre + stop_umber
stops.sort()



Frame_4 = collections.OrderedDict()
Frame_5 = collections.OrderedDict()
Frame_6 = collections.OrderedDict()

Seen_Stops = []
counter = 0
Stop_Stop = collections.OrderedDict()

Lengths = []

for stop in stops:
    Found = False
    #if stop not in Seen_Stops:
    Seen_Stops.append(stop)
    if stop != stops[-1]:
        for next_stop in stops[counter+1:]:
            diff = abs(stop-next_stop)
            Lengths.append(diff)

            if diff % 3 == 0:
                if next_stop not in Seen_Stops:
                    #if diff > 500:
                    frame = stop % 3
                    if frame == 0:
                        Frame_4.update({stop:next_stop})
                    if frame == 1:
                        Frame_5.update({stop:next_stop})
                    if frame == 2:
                        Frame_6.update({stop:next_stop})

                    Found = True
                    Seen_Stops.append(next_stop)
                    break

    # if Found == Fals:
    #     Stop_Stop

    counter +=1

print max(Lengths)
############################
Found_ORFs.write("Frame 4 \n")
for k, v in Frame_4.iteritems():
    Found_ORFs.write(str(k) + "\t" + "-" + "\t" + str(v) + "\n")
Found_ORFs.write("Frame 5 \n")
for k, v in Frame_5.iteritems():
    Found_ORFs.write(str(k) + "\t" + "-" + "\t" + str(v) + "\n")
Found_ORFs.write("Frame 6 \n")
for k, v in Frame_6.iteritems():
    Found_ORFs.write(str(k) + "\t" + "-" + "\t" + str(v) + "\n")
#############################
Longer_ORFs.write("Frame 4 \n")
for k, v in Frame_4.iteritems():
    if abs(k - v) > 200:
        Longer_ORFs.write(str(k) + "\t" + "-" + "\t" + str(v) + "\n")
Longer_ORFs.write("Frame 5 \n")
for k, v in Frame_5.iteritems():
    if abs(k - v) > 200:
        Longer_ORFs.write(str(k) + "\t" + "-" + "\t" + str(v) + "\n")
Longer_ORFs.write("Frame 6 \n")
for k, v in Frame_6.iteritems():
    if abs(k - v) > 200:
        Longer_ORFs.write(str(k) + "\t" + "-" + "\t" + str(v) + "\n")
##################################
Larger_ORFs.write("Frame 4 \n")
for k, v in Frame_4.iteritems():
    if abs(k - v) >= 350:
        Larger_ORFs.write(str(k) + "\t" + "-" + "\t" + str(v) + "\n")
Larger_ORFs.write("Frame 5 \n")
for k, v in Frame_5.iteritems():
    if abs(k - v) >= 350:
        Larger_ORFs.write(str(k) + "\t" + "-" + "\t" + str(v) + "\n")
Larger_ORFs.write("Frame 6 \n")
for k, v in Frame_6.iteritems():
    if abs(k - v) >= 350:
        Larger_ORFs.write(str(k) + "\t" + "-" + "\t" + str(v) + "\n")

