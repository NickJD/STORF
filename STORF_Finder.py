
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import numpy as np

import collections
import pdb

from DNA_Reverse_Complement import revCompIterative
import bisect
import sys

#inSequence = sys.argv[1]


try:

   #sys.argv[1:]
   minOrf = int(sys.argv[1]) #// do sth with sys,argv[1:]
except IndexError:
   minOrf = 10




sequence = ('Test0.txt')
#sequence = open('seq.txt',mode='rb')

Found_ORFs = open('./Predicted_ORFs.txt',mode='wb')
#Filtered_ORFs = open('./Filtered_Predicted_ORFs.txt',mode='wb')

#Matched_ORFs = open('./Matched_ORFs.txt',mode='wb')

#############################################################################

#pdb.set_trace()
orff_Positions = []
orff_Positions_Reverse = []

starts_stops = {}
starts_stops_within = []


def find_orfs_Forward(sequence):
    """ Finds all valid open reading frames in the string 'sequence', and
        returns them as a list"""

    stop_amber = find_all(sequence, 'TAG')
    stop_ochre = find_all(sequence, 'TAA')
    stop_umber = find_all(sequence, 'TGA')
    stops = stop_amber + stop_ochre + stop_umber
    stops.sort()


    ORF_Temps = collections.OrderedDict()


    seq_size = len(sequence)
    #
    # for f_stops in stops:
    #     f_t = seq_size - f_stops
    #     s = sequence[f_stops]
    #     f = f_t/3

    for stops_M in stops:
        Starts = []
        for st in stops:
            if (stops_M - st) % 3 == 0:
                Starts.append(st)
        if len(Starts) == 0 and stops[stops_M] not in Starts:
            Starts.append(stops[stops_M])


        p = min(Starts) +3
        #pos = ('1-' + str(p))
        seq = sequence[0:min(Starts) + 3]
        #pdb.set_trace()
        if len(seq) >= minOrf:
            ORF_Temps.update({1: p})



    for stops_M in reversed(stops):
        Ends = []
        for end in reversed(stops):
            if  (stops_M - end) % 3 == 0:
                Ends.append(end)
        if len(Ends) == 0:
            Ends.append(stops[-1])
        p = max(Ends)
        #pos = (str(p) + '-' + str(seq_size))
        seq = sequence[max(Ends):seq_size + 3]
        if len(seq) >= minOrf:
            ORF_Temps.update({p:seq_size})



    for stop_A in stops:
        for stop_B in stops:
            if stop_B > stop_A+10 and (stop_A - stop_B) % 3 == 0:  # Stop is in-frame and gap is more than 10 NT's

                #pos = (str(stop_A) + '-' + str(stop_B))
                seq = sequence[stop_A:stop_B + 3]
                if len(seq) >= minOrf:
                    ORF_Temps.update({stop_A: stop_B})
                break



    Frames = [0, 1, 2]  ### Sort frames out - Do I want 0-STOP or inframe maybe 1 - STOP??/
    stop_codons = ['TAG', 'TAA', 'TGA']
    for frame in Frames:
        t_seq = sequence[frame:]
        for i in range(3, len(t_seq) + 1, 3):
            codon = t_seq[i - 3:i]
            if codon in stop_codons:
                print i
                frame_seq = sequence[0:i + frame]
                if frame_seq not in ORF_Temps.itervalues():
                    if len(frame_seq) >= minOrf:
                        ORF_Temps.update({0:i})
                break



    Filtered_ORF = collections.OrderedDict()
    ORF_Temps = collections.OrderedDict(sorted(ORF_Temps.items()))
    for Start,Stop in ORF_Temps.items():

        if len(Filtered_ORF) == 0:
            Filtered_ORF.update({Start: Stop})
            continue

        F_Start = Filtered_ORF.keys()[-1]
        F_Stop = Filtered_ORF.get(F_Start)

        if F_Start == Start and F_Stop == Stop:
           break

        elif Stop == F_Stop and Start < F_Start:
            #if
            del Filtered_ORF[F_Start]
            Filtered_ORF.update({Start:Stop})

        elif Stop == F_Stop and Start > F_Start:
            continue

        elif Start >= F_Start and Stop <= F_Stop:
            continue

        elif Start == F_Stop:
            Filtered_ORF.update({Start: Stop})

        elif Start >= F_Start and Start < F_Stop:
            if abs(Start-F_Stop) <=10:
                Filtered_ORF.update({Start: Stop})
            elif abs(F_Start-Start) <=10:
                S = abs(Start-Stop)
                F = abs(F_Start-F_Stop)
                if S >= F:
                    del Filtered_ORF[F_Start]
                    Filtered_ORF.update({Start: Stop})
            elif Stop > F_Stop:

                Filtered_ORF.update({F_Start:Stop}) #Combine the two ORF's





        #elif {k:v for (k,v) in Filtered_ORF.items() if Start > v}:
        elif all(i < Start for i in Filtered_ORF.itervalues()):
            Filtered_ORF.update({Start:Stop})

            #
            # elif Stop == F_Stop and Start > F_Start:
            #     Covered = True
            #     break
            #
            # elif Start <= F_Start and Stop >= F_Stop:
            #     del Filtered_ORF[F_Start]
            #     Filtered_ORF.update({Start:Stop})
            #     Covered = True
            #
            # elif {k:v for (k,v) in Filtered_ORF.items() if Start > v}:
            #     Filtered_ORF.update({Start:Stop})
            #     Covered = True


        # if Covered == False:
        #
        #     if Start-100 <= F_Stop and Start >= F_Start and Stop > F_Stop+100 :
        #         #del Filtered_ORF[F_Pos]
        #         Filtered_ORF.update({Start:Stop})

    Final_ORFs = collections.OrderedDict()
    for F_Start, F_Stop in Filtered_ORF.items():
        pos = (str(F_Start) + '-' + str(F_Stop))
        seq = sequence[F_Start:F_Stop + 3]
        Final_ORFs.update({pos:seq})


    return Final_ORFs


def find_orfs_Reverse(sequence,selectedORFs):
    """ Finds all valid open reading frames in the string 'sequence', and
        returns them as a list"""

    stop_amber = find_all(sequence, 'TAG')
    stop_ochre = find_all(sequence, 'TAA')
    stop_umber = find_all(sequence, 'TGA')
    stops = stop_amber + stop_ochre + stop_umber
    stops.sort()


    ORF_Temps = collections.OrderedDict()


    seq_size = len(sequence)


    for stops_M in stops:
        Starts = []
        for st in stops:
            if (stops_M - st) % 3 == 0:
                Starts.append(st)
        if len(Starts) == 0 and stops[stops_M] not in Starts:
            Starts.append(stops[stops_M])

        p = min(Starts) +3
        pos = ('1-' + str(p))
        seq = sequence[0:min(Starts) + 3]
        if len(seq) >= minOrf:
            ORF_Temps.update({1: p})
    #break


    for stops_M in reversed(stops):
        Ends = []
        for end in reversed(stops):
            if  (stops_M - end) % 3 == 0:
                Ends.append(end)
        if len(Ends) == 0:
            Ends.append(stops[-1])
        p = min(Ends)
        pos = (str(p) + '-' + str(seq_size))
        seq = sequence[max(Ends):seq_size + 3]
        if len(seq) >= minOrf:
            ORF_Temps.update({p: seq_size})




    for stop_A in stops:
        for stop_B in stops:
            if stop_B > stop_A + 10 and (stop_A - stop_B) % 3 == 0:  # Stop is in-frame and gap is more than 10 NT's

                pos = (str(stop_A) + '-' + str(stop_B))
                seq = sequence[stop_A:stop_B + 3]
                if len(seq) >= minOrf:
                    ORF_Temps.update({stop_A: stop_B})
                break






    Frames = [0, 1, 2]  ### Sort frames out - Do I want 0-STOP or inframe maybe 1 - STOP??/
    stop_codons = ['TAG', 'TAA', 'TGA']
    for frame in Frames:
        t_seq = sequence[frame:]
        for i in range(3, len(t_seq) + 1, 3):
            codon = t_seq[i - 3:i]
            if codon in stop_codons:
                print i
                frame_seq = sequence[0:i + frame]
                if frame_seq not in ORF_Temps.itervalues():
                    if len(frame_seq) >= minOrf:
                        ORF_Temps.update({0:i})

                break

                         ###Check that first direction works.

    if FilterReverse == True: #Filter with Anti-Sense/Sense
        readLen = len(sequence)#Sequence Length
        Final_Orfs = collections.OrderedDict()
        for r_start, r_stop in ORF_Temps.iteritems(): #Get Anti-Sense ORF Positions
            # r_start_a = k_pos.split("-")[0]
            # r_stop_a = k_pos.split("-")[1]


            #pos = str(r_start_a) + '-' + str(r_stop_a)
            #r_stop = int(r_stop_a)-3

            seq = sequence[int(r_start):int(r_stop) + 3]
            if len(seq) >=10:
                Final_Orfs.update({r_start:r_stop})


    Filtered_ORF = collections.OrderedDict()
    for pos in selectedORFs.keys():
        Start = int(pos.split('-')[0])
        Stop = int(pos.split('-')[1])

        Final_Orfs.update({Start:Stop})
    Final_Orfs = collections.OrderedDict(sorted(Final_Orfs.items()))
    for Start,Stop in Final_Orfs.items():

        #Start = int(pos.split('-')[0])
        #Stop = int(pos.split('-')[1])
        if len(Filtered_ORF) == 0:
            Filtered_ORF.update({Start: Stop})
            continue


            # pos = (str(Start) + '-' + str(Stop))

        for F_Start, F_Stop in Filtered_ORF.items():

            if F_Start == Start and F_Stop == Stop:
                break

            elif Stop == F_Stop and Start < F_Start:
                del Filtered_ORF[F_Start]
                Filtered_ORF.update({Start: Stop})

            # elif {k:v for (k,v) in Filtered_ORF.items() if Start > v}:
            elif all(i < Start for i in Filtered_ORF.itervalues()):
                Filtered_ORF.update({Start: Stop})

                #
                # elif Stop == F_Stop and Start > F_Start:
                #     Covered = True
                #     break
                #
                # elif Start <= F_Start and Stop >= F_Stop:
                #     del Filtered_ORF[F_Start]
                #     Filtered_ORF.update({Start:Stop})
                #     Covered = True
                #
                # elif {k:v for (k,v) in Filtered_ORF.items() if Start > v}:
                #     Filtered_ORF.update({Start:Stop})
                #     Covered = True


                # if Covered == False:
                #
                #     if Start-100 <= F_Stop and Start >= F_Start and Stop > F_Stop+100 :
                #         #del Filtered_ORF[F_Pos]
                #         Filtered_ORF.update({Start:Stop})

    Final_ORFs = collections.OrderedDict()
    for F_Start, F_Stop in Filtered_ORF.items():
        pos = (str(F_Start) + '-' + str(F_Stop))
        seq = sequence[F_Start:F_Stop + 3]
        Final_ORFs.update({pos: seq})

        return Final_ORFs

        # pos = (str(Start) + '-' + str(Stop))
        # for F_Pos, F_Seq in Filtered_ORF.items():
        #     F_Start = int(F_Pos.split('-')[0])
        #     F_Stop = int(F_Pos.split('-')[1])
        #
        #     if Stop == F_Stop and Start < F_Start:
        #         del Filtered_ORF[F_Pos]
        #         Filtered_ORF.update({pos:seq})
        #
        #     elif F_Start == Start and F_Stop == Stop:
        #         continue
        #
        #     elif Start <= F_Start and Stop >= F_Stop:
        #         del Filtered_ORF[F_Pos]
        #         Filtered_ORF.update({pos: seq})
        #
        #     elif Start-100 <= F_Stop and Start >= F_Start and Stop > F_Stop :
        #         #del Filtered_ORF[F_Pos]
        #         Filtered_ORF.update({pos: seq})

    #return ORF_Temps
    #return Filtered_ORF

    #return Final_Orfs






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
with open(sequence) as f:
    lines = f.readlines()

    for line in lines:
        if line is lines[-1]:
            Contigs.append(Contig)
            Contig = ""

        line = line.replace("\n", "")
        if ">" not in line:
            Contig += str(line)
        elif ">" in line:
            Contigs.append(Contig)
            Contig = ""
        elif line is None:
            print "end"
            Contigs.append(Contig)
            Contig = ""

#genedict = []
#genedict = seq.replace("\n","")





orfdict = {}

all_orfs = []
all_selectedORFs = []
all_orfs_rev = []
all_selectedORFs_r = []





FilterReverse = False

selectedORFs_f = collections.OrderedDict()
selectedORFs_r = collections.OrderedDict()

Conum = 1

Contigs = filter(None, Contigs)

for con in Contigs:
    if con != '':

        temp_ORFs_f =  find_orfs_Forward(con) #Get Forward ORFs
        selectedORFs_f = temp_ORFs_f



        con_rev = revCompIterative(con) #Get Reverse Compliment

        FilterReverse = True

        temp_ORFs_r = find_orfs_Reverse(con_rev,selectedORFs_f) # Get Reverse ORFs
        selectedORFs_r.update(temp_ORFs_r)






###############################################
    Found_ORFs.write("ORFs for Contig " + str(Conum) + "\n")
    Found_ORFs.write("Forward Strand \n")
    for k_pos, v_seq in selectedORFs_f.iteritems():
        r_start = k_pos.split("-")[0]
        r_stop = k_pos.split("-")[1]
        Found_ORFs.write(str(k_pos) + "\t" + "+" + "\t" + v_seq + "\n")


    Found_ORFs.write("Revserse Compliment \n")
    for k_pos, v_seq in selectedORFs_r.iteritems():
        r_start = k_pos.split("-")[0]
        r_stop = k_pos.split("-")[1]
        Found_ORFs.write(str(k_pos) + "\t" + "-" + "\t" + v_seq + "\n")




