import sys
import numpy as np


sequence = sys.argv[1]
minOrfSize = int(sys.argv[2])
maxOrfSize = int(sys.argv[3])
stop_Codons = sys.argv[4]
output = sys.argv[5]
import collections


def revCompIterative(watson): #Gets Reverse Complement
    complements = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C', 'N': 'N'}
    watson = watson.upper()
    watsonrev = watson[::-1]

    crick = ""


    for nt in watsonrev:
        crick += complements[nt]

    return crick


def find_all(sequence, subsequence): # Finds all Stop locations
    ''' Returns a list of indexes within sequence that are the start of subsequence'''
    start = 0
    idxs = []
    next_idx = sequence.find(subsequence, start)

    while next_idx != -1:
        idxs.append(next_idx)
        start = next_idx + 1  # Move past this on the next time around
        next_idx = sequence.find(subsequence, start)

    return idxs


################### NT - AA


gencode = {
      'ATA':'I', 'ATC':'I', 'ATT':'I', 'ATG':'M',
      'ACA':'T', 'ACC':'T', 'ACG':'T', 'ACT':'T',
      'AAC':'N', 'AAT':'N', 'AAA':'K', 'AAG':'K',
      'AGC':'S', 'AGT':'S', 'AGA':'R', 'AGG':'R',
      'CTA':'L', 'CTC':'L', 'CTG':'L', 'CTT':'L',
      'CCA':'P', 'CCC':'P', 'CCG':'P', 'CCT':'P',
      'CAC':'H', 'CAT':'H', 'CAA':'Q', 'CAG':'Q',
      'CGA':'R', 'CGC':'R', 'CGG':'R', 'CGT':'R',
      'GTA':'V', 'GTC':'V', 'GTG':'V', 'GTT':'V',
      'GCA':'A', 'GCC':'A', 'GCG':'A', 'GCT':'A',
      'GAC':'D', 'GAT':'D', 'GAA':'E', 'GAG':'E',
      'GGA':'G', 'GGC':'G', 'GGG':'G', 'GGT':'G',
      'TCA':'S', 'TCC':'S', 'TCG':'S', 'TCT':'S',
      'TTC':'F', 'TTT':'F', 'TTA':'L', 'TTG':'L',
      'TAC':'Y', 'TAT':'Y', 'TAA':'', 'TAG':'',
      'TGC':'C', 'TGT':'C', 'TGA':'', 'TGG':'W'}

def translate_frameshifted( sequence ):
      translate = ''.join([gencode.get(sequence[3*i:3*i+3],'X') for i in range(len(sequence)//3)])
      return translate



#############################










Genome = ''

Contigs = collections.OrderedDict()

with open (sequence, 'rb') as Gfile: #Input file
    lines = Gfile.read().splitlines()
    Contig = ""
    Contig_ID = ""
    for line in lines:
        line = line.replace("\n", "")
        if line is lines[-1]:
            Contigs.update({Contig_ID:Contig})
        if ">" in line and line is not lines[0]:
            Contigs.update({Contig_ID: Contig})
            Contig = ""
            Contig_ID = line
        elif ">" in line:
            Contig = ""
            Contig_ID = line

        if ">" not in line:
            Contig += str(line)


for Contig_ID, Contig in Contigs.items():



    stops = []
    for stop_codon in stop_Codons.split(','):
        stops.extend(find_all(Contig, stop_codon))


    stops.sort()

    Genome_Size = len(Contig)


    Frames = collections.OrderedDict()



    Seen_Stops = []
    counter = 0
    Stop_Stop = collections.OrderedDict()
    ORFNum = 0
    Lengths = []

    for stop in stops: # Finds Stop-Stop#
        GeneA = np.zeros((Genome_Size), dtype=np.int)
        Found = False
        Seen_Stops.append(stop)
        if stop != stops[-1]:
            for next_stop in stops[counter+1:]:
                diff = abs(stop-next_stop)
                Lengths.append(diff)

                if diff % 3 == 0:
                    if next_stop not in Seen_Stops:
                        if diff >= minOrfSize and diff <= maxOrfSize:

                            if ORFNum >1:
                                GeneA[stop:next_stop] = [1] * (next_stop - stop)
                                length = next_stop - stop
                                prevlength = PrevNext - PrevStop
                                Gene_AND = GenePrev & GeneA

                                if stop > PrevStop and next_stop < PrevNext:
                                    break

                                if np.count_nonzero(Gene_AND) <= (length * 0.5):



                                    Frames.update({stop:next_stop})
                                    Found = True
                                    Seen_Stops.append(next_stop)
                            else:
                                Frames.update({stop:next_stop})
                                Found = True
                                Seen_Stops.append(next_stop)
                        GenePrev = GeneA
                        PrevStop = stop
                        PrevNext = next_stop
                        ORFNum += 1
                        break

        counter +=1




    #################################

    with open(output,'a') as out:
        out.truncate()
        #out.write("Frame 1 \n")
        STORF = 0
        for k, v in Frames.iteritems():
            Contig_ID = Contig_ID.replace('>','')
            out.write(">"+str(STORF)+":"+str(Contig_ID)+','+str(k) + "\t" + "+" + "\t" + str(v) + "\n")
            #Amino = translate_frameshifted(Contig[k:v+3])
            out.write(Contig[k:v+3]+'\n')
            STORF +=1























    ###### Reversed

    Contig_rev = revCompIterative(Contig)



    stops = []
    for stop_codon in stop_Codons.split(','):
        stops.extend(find_all(Contig_rev, stop_codon))
    stops.sort()



    Frames_rev = collections.OrderedDict()


    Seen_Stops = []
    counter = 0
    Stop_Stop = collections.OrderedDict()

    Lengths = []

    for stop in stops:
        Found = False
        Seen_Stops.append(stop)
        if stop != stops[-1]:
            for next_stop in stops[counter+1:]:
                diff = abs(stop-next_stop)
                Lengths.append(diff)

                if diff % 3 == 0:
                    if next_stop not in Seen_Stops:
                        if diff >= minOrfSize and diff <= maxOrfSize:
                            Frames_rev.update({stop:next_stop})


                            Found = True
                            Seen_Stops.append(next_stop)
                        break



        counter +=1


    ############################




    ###############
    with open(output,'a') as out:
        #out.write("Frame 4 \n")
        STORF = 0
        for k, v in Frames_rev.iteritems():
            Contig_ID = Contig_ID.replace('>', '')
            out.write(">"+str(STORF)+":"+str(Contig_ID)+','+str(k) + "\t" + "-" + "\t" + str(v) + "\n")
            #Amino = translate_frameshifted(Contig_rev[k:v+3])
            out.write(Contig_rev[k:v+3]+'\n')
            STORF += 1
