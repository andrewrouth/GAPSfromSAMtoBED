import re
import argparse
parser = argparse.ArgumentParser()
parser.add_argument("SAMIN", help="Description")
parser.add_argument("Output", help="Description")
args = parser.parse_args()
SAMIN = str(args.SAMIN)
Outputfile = str(args.Output)

##      -----------------------------------------------------------------------------------------------------------------

def chunks(l, n):
     ##Yield successive n-sized chunks from l##
     for i in range(0, len(l), n):
         yield l[i:i+n]
         
class BED(object):    
    def __init__(self, Name, BEDEntry):
        self.Name = Name
        self.Ref = BEDEntry[0]
        self.From = BEDEntry[1]
        self.To = BEDEntry[2]
        self.Name = BEDEntry[3]
        self.Count = 1
        self.Dir = BEDEntry[5]

def FindDelCoords(StartNuc, MappingStats, Ref):
    ##StartNuc in 1-pos
    Countable = ['M', 'X', 'N', 'D']
    Events = []
    Events.append(str(StartNuc))
    MappingStats = re.split('(\D+)', MappingStats)[:-1]
    NLocs = [int(i/2) for i, x in enumerate(MappingStats) if x=='D' or x=='N']
    MappingStats = list(chunks(MappingStats,2))
    n=0
    ##Gaps
    if NLocs:
        for i in NLocs:
            Segment = MappingStats[0:i]
            SegmentLength = 0
            DelLength = int(MappingStats[i][0])
            #All_Lengths.append(DelLength)
            for j in Segment:
                if j[1] in Countable:
                    SegmentLength += int(j[0])
            FromCoord = StartNuc + SegmentLength
            ToCoord = FromCoord + DelLength - 1
            #Event = str(FromCoord) + '^' + str(ToCoord)
            BEDEntry = [Ref, str(FromCoord), str(ToCoord), 'Deletion', '1', Dir]
            BEDName = '\t'.join(BEDEntry)
            n+=1
            if BEDName not in BEDDict:
                BEDDict[BEDName] = BED(BEDName, BEDEntry)
            else:
                BEDDict[BEDName].Count += 1
    ##FindEnd
    SegmentLength = 0
    for j in MappingStats:
        if j[1] in Countable:
            SegmentLength += int(j[0])
    EndNuc = StartNuc + SegmentLength - 1
    Events.append(str(EndNuc))
    return Events, n

#exampleread = ['SNL122:189:HGHNTBCXY:1:1107:14825:52979',
#                 '0',
#                 'NC_004146.1_FHV_RNA1.seq',
#                 '1082',
#                 '255',
#                 '32M1X27M52N32M1521N52M',
#                 '*',
#                 '0',
#                 '0',
#                 'TATCCGCCACCCAATCTGTCAACGCTAGGCTTGTCGGTATGGGACACAAGGACCCGCAATTAGTCCAACTGTGTATAAACCTACAATGCCACAAACTCGCGCTAATCCAGGAACTTCCCGACCGCATTCAAACGGCGGTGGAAG',
#                 'IIIIIIIIIIIIHIIIIIIIIIIIIIHIIHIIIIIIIIIIHIIHHIIIHHHIIIIIIIIIIIIIIIHIHHHIIIIHIIGHIIIIHIIIIIIIIIIFHIIHDIIIGIIIIIIIGIIIIIIIIIIHIIIHIIIIIIIIHIIIGIII',
#                 'NM:i:1']
BEDDict = {}
with open(SAMIN) as SAMIn:
    line = SAMIn.readline()
    while line:
        if line[0] != '@':
            data = line.split()
            Name = data[0]
            flag = bin(int(data[1]))
            if data[1] != '4':
                StartNuc = int(data[3])
                Ref = data[2]
                if len(flag) > 6 and flag[-5] == '1':
                    #A flag of '16' means the read mapped to the -ve sense ref
                    Dir = '-'
                else:
                    Dir = '+'
                MappingStats = data[5]
                Events, n = [], 0
                Events, n = FindDelCoords(StartNuc, MappingStats, Ref)
#                if Name == exampleread[0]:
#                    print(Events, n)
            else:
                    pass
        else:
                pass
        line = SAMIn.readline()
        
#BEDFILE FORMAT:
#track name=Virus_Recombinations description="Virus_Recombinations" graphType=junctions
#REF	From	To 	Name	Count	Dir
#e.g.
#NC_004146.1_FHV_RNA1.seq	1242	2326	Deletion	486	+

FinalBEDList = []
for j in BEDDict:
    FinalBEDList.append([BEDDict[j].Ref, BEDDict[j].From, BEDDict[j].To, BEDDict[j].Name, str(BEDDict[j].Count), BEDDict[j].Dir])
FinalBEDList.sort(key=lambda a:int(a[4]), reverse=True)

with open(Outputfile,'w') as Out:
    Out.write('track name=Virus_Recombinations description="Virus_Recombinations" graphType=junctions\n')
    for j in FinalBEDList:
        Out.write('\t'.join(j) + '\n')