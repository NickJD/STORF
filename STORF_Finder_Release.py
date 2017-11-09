import sys
#sequence = open('seq.txt',mode='rb')
from DNA_Reverse_Compliment import revCompIterative


sequence = sys.argv[1]
orfSize = int(sys.argv[2])
stop_Codons = sys.argv[3]
output = sys.argv[4]
import re
import collections



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




Genome = ''


with open (sequence, 'rb') as Gfile:
    for line in Gfile:
        line = line.replace("\n","")
        if ">" not in line:
            Genome += str(line)

stops = []
for stop_codon in stop_Codons.split(','):
    stops.extend(find_all(Genome, stop_codon))


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

                    if diff >= orfSize:
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
with open(output,'w') as out:
    out.truncate()
    out.write("Frame 1 \n")
    for k, v in Frame_1.iteritems():
        out.write(str(k) + "\t" + "+" + "\t" + str(v) + "\n")
    out.write("Frame 2 \n")
    for k, v in Frame_2.iteritems():
        out.write(str(k) + "\t" + "+" + "\t" + str(v) + "\n")
    out.write("Frame 3 \n")
    for k, v in Frame_3.iteritems():
        out.write(str(k) + "\t" + "+" + "\t" + str(v) + "\n")

###### Reversed
Genome_rev = revCompIterative(Genome)



stops = []
for stop_codon in stop_Codons.split(','):
    stops.extend(find_all(Genome_rev, stop_codon))
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
                    if diff >= orfSize:
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
with open(output,'a') as out:
    out.write("Frame 4 \n")
    for k, v in Frame_4.iteritems():
        out.write(str(k) + "\t" + "-" + "\t" + str(v) + "\n")
    out.write("Frame 5 \n")
    for k, v in Frame_5.iteritems():
        out.write(str(k) + "\t" + "-" + "\t" + str(v) + "\n")
    out.write("Frame 6 \n")
    for k, v in Frame_6.iteritems():
        out.write(str(k) + "\t" + "-" + "\t" + str(v) + "\n")