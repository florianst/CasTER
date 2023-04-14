import os
import re
import subprocess
import pickle
import numpy as np
import pandas as pd
import time
import threading
import queue
from shutil import copyfile
from tqdm import tqdm

class Oligo:
    def __init__(self, oligoNumber, guide_sequence, target_sequence, pam_index, strand, totalReads):
        self.oligoNumber, self.totalReads = oligoNumber, totalReads
        self.guide_sequence, self.target_sequence, self.pam_index, self.strand = guide_sequence, target_sequence, pam_index, strand

        self.dict_outcomes = {}

    def add_outcome(self, isInsertion, indelLength, indelStart, indelStartRel, relFreq, outcome_sequence):
        # TODO: dictionary keys should be (indelLength, indelStart, isInsertion)!!
        # otherwise we overwrite (indelLength, indelStart) if it appears for both insertion and deletion
        self.dict_outcomes[(indelLength, indelStart)] = [isInsertion, indelStartRel, relFreq, outcome_sequence, None]

    def add_energies(self, potIndelLength, potIndelStart, mosaics_energies):
        # fill in relFreq appropriately
        if (potIndelLength, potIndelStart) in self.dict_outcomes.keys():
            self.dict_outcomes[(potIndelLength, potIndelStart)][-1] = mosaics_energies

        else:
            isInsertion = False
            potIndelStartRel = potIndelStart - self.pam_index
            relFreq, outcome_sequence = 0, '' # fill relFreq with 0 for unmeasured outcomes
            self.dict_outcomes[(potIndelLength, potIndelStart)] = [isInsertion, potIndelStartRel, relFreq, outcome_sequence, mosaics_energies]

    def get_df(self):
        # sanity checks on dict_outcomes
        if len(self.dict_outcomes) < 1: return None
        for k, v in self.dict_outcomes.items():
            if v[-2] is None: # deal with a legacy version of the Oligo class
                self.dict_outcomes[k] = self.dict_outcomes[k][:-2]+[self.dict_outcomes[k][-1]]

        df = pd.DataFrame(data=self.dict_outcomes).transpose().reset_index()
        df.columns=["indelLength", "indelStart", "isInsertion", "potIndelStartRel", "relFreq", "outcome_sequence", "mosaics_energies"]
        df = pd.concat([df.drop(['mosaics_energies'], axis=1), df['mosaics_energies'].apply(pd.Series)], axis=1)

        energy_columns = [c for c in df.columns if "energy" in c]
        for c in energy_columns:
            df[c] = df[c].apply(lambda d: d if isinstance(d, list)       else [np.nan, None, np.nan])
            df = df.join(pd.DataFrame(df[c].tolist(), index=df.index, columns=["val", "seq", "bppt"]).add_prefix(c+'_'))
            df.drop(c, axis=1, inplace=True)
            
        return df
    

# KL Divergence between the two indel profiles (non-symmetric) - from felicityallen/SelfTarget
def KL(p1, p2, ignore_null=True, missing_count=0.5):
    
    p1_indels = set([x for x in p1 if p1[x]>0]) #and (x != '-' or not ignore_null)
    p2_indels = set([x for x in p2 if p2[x]>0])
    common = p1_indels.intersection(p2_indels)
    p1_only = p1_indels.difference(p2_indels)
    p2_only = p2_indels.difference(p1_indels)
    
    p1_total = sum([p1[x] for x in p1_indels]) + missing_count*len(p2_only)
    p2_total = sum([p2[x] for x in p2_indels]) + missing_count*len(p1_only)
    
    if p1_total > 0 and p2_total > 0:
        norm1, norm2 = 1.0/p1_total, 1.0/p2_total
        score = 0.0
        for indel in common:
            score += p1[indel]*norm1*np.log2(p1[indel]*norm1/(p2[indel]*norm2))
        for indel in p1_only:
            score += p1[indel]*norm1*np.log2(p1[indel]*norm1/(missing_count*norm2))
        for indel in p2_only:
            score += missing_count*norm1*np.log2(missing_count*norm1/(p2[indel]*norm2))
    else: score = np.nan
    return score

def symmetricKL(profile1, profile2, ignore_null=True):
    return 0.5*KL(profile1, profile2, ignore_null) + 0.5*KL(profile2, profile1, ignore_null)

def averageKL(df_1, df_2):
    df_1 = df_1.copy() * 1000 # achieve comparability with FORECasT - they normalise total read counts to 1000
    df_2 = df_2.copy() * 1000
    dict_1, dict_2 = df_1.to_dict(orient='records'), df_2.to_dict(orient='records')
    symmetric_KL_dist = [symmetricKL(p1, p2) for p1, p2 in zip(dict_1, dict_2)]

    return np.median(symmetric_KL_dist), symmetric_KL_dist


def changeLine(linenumber, file, newline): # change single line in file
    with open(file, 'r') as f:
        lines = f.readlines()

    with open(file, 'w') as f:
        for i, line in enumerate(lines):
            if (i == linenumber): line = newline+'\n'
            f.write(line)

DNAbaseCorresp = {ord('A'):'T', ord('T'):'A', ord('C'):'G', ord('G'):'C'}

def getOppositeBases(sequence): # returns opposite base pairs of given sequence
    sequence = sequence.upper().translate(DNAbaseCorresp)
    return sequence

def parseIndelIdentifier(identifier):
    identifierLength, identifierPos = identifier.split('_')[:2]
    identifierLength, identifierPos = re.split('([+-]?\d+)', identifierLength), re.split('([+-]?\d+)', identifierPos) # separate signed integers and letters
    if (identifierLength[0] == "I"): isInsertion = True
    else:                            isInsertion = False
    indelLength    = int(identifierLength[1]) # length of indel
    indelStartRel  = int(identifierPos[1])-2  # when subtracting two, this becomes the position of the last unaltered nucleotide relative to the beginning of the PAM
    
    return isInsertion, indelLength, indelStartRel

def find_microhomologies(left, right):
    mhs, offsets = [], []
    # shift right part underneath left
    for right_offset in np.arange(-(len(right)-1), len(left)):
        mh = []
        for idx in range(min(len(left), len(right)-right_offset)):
            if left[idx] == right[max(0, right_offset+idx)]:
                mh.append(idx)
            else:
                if len(mh) > 1: 
                    offsets.append(-right_offset)
                    mhs.append(mh)
                mh = []
        if len(mh) > 1: 
            offsets.append(-right_offset)
            mhs.append(mh)

    return zip(mhs, offsets)


def loadOligos(lookupFile): # open oligo information from the lookup file (S3) given in Allen et al
    oligos = {}
    with open(lookupFile) as fp:
        for line in fp:
            if (line[:5] == "Oligo"): 
                line = line.strip().split('\t')
                oligoNumber = int(re.sub("\D", "", line[0])) # strip everything apart from numbers
                guide_sequence, target_sequence, pam_index, strand = line[1], line[2], int(line[6]), line[7] == "FORWARD"
                if (strand == False): target_sequence = getOppositeBases(target_sequence)[1:] # need to invert bases in case strand is reverse
                
                oligos[oligoNumber] = Oligo(oligoNumber, guide_sequence, target_sequence, pam_index, strand, -1) # totalReads is filled in later
    
    return oligos


def loadAllenData(oligoDir, lookupFile, maxOligos=-1, mosaicsWindow=21, numWorkers=8, maxIndelLength=20):
    # read lookup file for guide and target sequences
    oligos = loadOligos(lookupFile)
    
    global data, mosaics_errorcount, progress_bar, indelLengthRange, mosaics_energy_cache

    indelLengthRange = maxIndelLength # calculate deletions of length [0, indelLengthRange]

    mosaics_errorcount = 0
    data = []
    if mosaics_energy_cache is None: mosaics_energy_cache = {}
    q = queue.Queue()
    currentOligoNum = -1
    for (dirpath, _dirnames, filenames) in os.walk(oligoDir):
        for file in filenames:
            with open(os.path.join(dirpath, file)) as fp:
                for line in fp:
                    line = line.strip('\n')
                    if (line[:3] == '@@@'):
                        if (currentOligoNum > -1): oligos[currentOligoNum].totalReads = totalReads
                        totalReads = 0
                        currentOligoNum = int(re.sub("\D", "", line.split('@')[-1]))
                        # sum observations for each oligo and write observation frequencies into oligos dictionary retroactively
                        continue
                    
                    line = line.split('\t')
                    if (len(line) < 3): continue
                    
                    totalReads += int(line[1])
                else: 
                    oligos[currentOligoNum].totalReads = totalReads # makes sure last oligo in each file is saved as well
                    totalReads = 0
            
            # go through same file again, this time filling in all outcomes for a given oligo
            with open(os.path.join(dirpath, file)) as fp:
                for line in fp:
                    line = line.strip('\n')
                    if (line[:3] == '@@@'): 
                        currentOligoNum = int(re.sub("\D", "", line.split('@')[-1]))
                        continue
    
                    line = line.split('\t')
                    if (len(line) < 3): continue
                    
                    identifier, mutation_count, outcome_sequence = line[:3]
                    
                    isInsertion, indelLength, indelStartRel = parseIndelIdentifier(identifier)
                    indelStart = indelStartRel + oligos[currentOligoNum].pam_index
                    
                    relFreq = float(mutation_count)/oligos[currentOligoNum].totalReads

                    oligos[currentOligoNum].add_outcome(isInsertion, indelLength, indelStart, indelStartRel, relFreq, outcome_sequence)
            
                if (os.path.exists('atompos.out')): os.remove('atompos.out') # remove file where final atom positions are saved - tends to get quite big, also jumbled by the workers

    # prepare MOSAICS process queue: one process per oligo
    oligo_objects = list(oligos.values())
    np.random.shuffle(oligo_objects)
    for currentOligo in oligo_objects:
        q.put([currentOligo, mosaicsWindow]) # if len(currentOligo.dict_outcomes) > 0: 
        if q.qsize() >= maxOligos and maxOligos > 0: break
    
    # start threads for a fixed number of workers, they write their results in data array
    progress_bar = tqdm(total=q.qsize()*(indelLengthRange+2)*(indelLengthRange+1)//2)
    threads = [threading.Thread(target=MOSAICSworker, args=(q,)) for _i in range(numWorkers)]
    for thread in threads:
        thread.start() # each thread takes one array from the Queue, executes MOSAICS and writes result into data array, then takes next array
        q.put(None)  # one EOF marker for each thread
        
    for thread in threads:
        thread.join() # wait until the workers have processed everything in the queue
    progress_bar.close()

    return data, mosaics_errorcount                 


def MOSAICSworker(q):
    global data, mosaics_errorcount, progress_bar, indelLengthRange, mosaics_energy_cache
    padding = 20 # number of bases padding around cut site

    while True:
        # TODO: measure execution time, try to minimise

        params = q.get()
        if params is None: return # end worker at last element in queue (None)
        
        # loop through all outcomes in the currentOligo, retrieve all relevant MOSAICS energies, save as single training point
        currentOligo, mosaicsWindow = params

        cutsite = currentOligo.pam_index-3

        temp_path = '' #'/home/ssd/'
        # write input file with specific files so that multiple workers can work in same directory
        inputFile = temp_path+'temp/input_'+str(currentOligo.oligoNumber)
        pdbFile   = temp_path+'temp/run_'+str(currentOligo.oligoNumber)+'.pdb'
        tempFile  = temp_path+'temp/temp_'+str(currentOligo.oligoNumber)+'.pdb'
        epotFile  = temp_path+'temp/epot_'+str(currentOligo.oligoNumber)+'.out'
        copyfile('input', inputFile)
        changeLine(49, inputFile, '  \epot_file{'+epotFile+'}')
        changeLine(46, inputFile, '  \pos_init_file{'+pdbFile+'}')

        # enumerate all deletions and start positions
        for potIndelLength in range(indelLengthRange+1):
            # augment Allen data with MOSAICS conjugate gradient energy of not yet unwound system
            intervalLength = mosaicsWindow // 2
            
            indelLengthAdjusted = potIndelLength
            while (indelLengthAdjusted >  10.5): indelLengthAdjusted -= 10.5 # this is basically like modulo but not quite
            while (indelLengthAdjusted < -10.5): indelLengthAdjusted += 10.5
            leftInTwoTurns = 2*10.5-indelLengthAdjusted
            leftInMosaicsWindow = leftInTwoTurns*(1-abs((mosaicsWindow/2-10.5))/(mosaicsWindow/2)) # correction factor such that the MOSAICS window takes up the torsion
            
            while (leftInMosaicsWindow > 12.5): leftInMosaicsWindow /= 2  # 12.5 bppt is an arbitrary value here (the highest the current pipeline can reliably produce the correct atom numbers with), should be more like 15
            bppt = leftInMosaicsWindow
            #if (leftInMosaicsWindow == 12.5): bppts = [leftInMosaicsWindow/2, leftInMosaicsWindow]

            isInsertion = False # TODO: deal with insertions

            for potIndelStart in np.arange(max(max(0, cutsite-potIndelLength-1), cutsite-padding), min(min(cutsite, len(currentOligo.target_sequence)-potIndelLength), cutsite+padding-potIndelLength)):
                # is the current outcome caused by a microhomology?
                # knowing mh_length, we can get the MMEJ sequences
                for i in list(range(10))[::-1]:
                    mh_length = i
                    if (currentOligo.target_sequence[potIndelStart-i:potIndelStart] == currentOligo.target_sequence[potIndelStart+potIndelLength:potIndelStart+potIndelLength+i]): break

                mosaicsSequence = currentOligo.target_sequence[potIndelStart+potIndelLength//2-intervalLength:potIndelStart+potIndelLength//2+intervalLength+1]
                if (len(mosaicsSequence) < intervalLength): 
                    print(mutation_id, bppt, mosaicsSequence, "sequence too short")
                    mosaics_errorcount += 1
                    continue # ignore sequences which do not fill the whole intervalLength (oligo target_sequence too short)

                sequence_X = currentOligo.target_sequence[max(0, cutsite-padding):potIndelStart]
                sequence_Y = currentOligo.target_sequence[potIndelStart:potIndelStart+potIndelLength]
                sequence_Z = currentOligo.target_sequence[potIndelStart+potIndelLength:min(cutsite+padding, len(currentOligo.target_sequence))]

                # torsion energies
                mosaics_runs = [['mosaics_energy',           mosaicsSequence,       bppt, 0, 0], 
                                ['mosaics_energy_noTorsion', mosaicsSequence,       10.5, 0, 0],
                                ['mosaics_energy_XZ',        sequence_X+sequence_Z, 10.5, 0, 0],
                                ['mosaics_energy_Y',         sequence_Y,            10.5, 0, 0]
                                ]

                # MMEJ energies
                if mh_length > 1 and not isInsertion: 
                    pad = mh_length if potIndelStart > cutsite else 0
                    repaired_sequence = currentOligo.target_sequence[cutsite-padding:min(potIndelStart, cutsite)] + currentOligo.target_sequence[potIndelStart+potIndelLength:cutsite+padding]
                    mosaics_runs.extend([
                        ['mosaics_energy_R',  currentOligo.target_sequence[cutsite-padding          :cutsite],                                10.5, 0,                                  0],
                        ['mosaics_energy_L',  currentOligo.target_sequence[cutsite                  :cutsite+padding],                        10.5, 0,                                  0],
                        ['mosaics_energy_R2', currentOligo.target_sequence[cutsite-padding          :min(potIndelStart, cutsite)],            10.5, 0,                                  mh_length],
                        ['mosaics_energy_r2', currentOligo.target_sequence[potIndelStart-mh_length  :cutsite],                                10.5, cutsite-(potIndelStart-mh_length),  0],
                        ['mosaics_energy_r1', currentOligo.target_sequence[potIndelStart            :cutsite],                                10.5, 0,                                  cutsite-potIndelStart],
                        ['mosaics_energy_L2', currentOligo.target_sequence[potIndelStart+potIndelLength:cutsite+padding],                     10.5, mh_length,                          0],
                        ['mosaics_energy_l2', currentOligo.target_sequence[cutsite                  :potIndelStart+potIndelLength+pad],       10.5, potIndelStart+potIndelLength+pad-cutsite, 0],
                        ['mosaics_energy_l1', currentOligo.target_sequence[cutsite+pad              :potIndelStart+potIndelLength+mh_length], 10.5, 0,                                  potIndelStart+potIndelLength+mh_length-(cutsite+pad)],
                        ['mosaics_energy_RL2',repaired_sequence,                                                                              10.5, 0,                                  0],
                    ])

                open(epotFile, 'a').close() # create epot file

                mosaics_energies, runs_success = {}, True
                for run_spec in mosaics_runs:
                    run_name, run_sequence, run_bppt, remove_bases_A, remove_bases_B = run_spec

                    if run_sequence != '':
                        # check if energy is in cache
                        if (run_sequence, run_bppt, remove_bases_A, remove_bases_B) in mosaics_energy_cache.keys():
                            mosaics_energies[run_name] = [mosaics_energy_cache[(run_sequence, run_bppt, remove_bases_A, remove_bases_B)], run_sequence, run_bppt]
                            continue

                        # if not, run MOSAICS
                        process = subprocess.Popen('./pipeline.sh '+run_sequence+' '+str(run_bppt)+' '+str(remove_bases_A)+' '+str(remove_bases_B)+' '+inputFile+' '+pdbFile+' '+tempFile, shell=True, stdout=subprocess.PIPE)
                        process.wait()
                        subprocess_output = str(process.communicate([0])[0], 'utf-8').split('\n')[-3]
                    
                        # read energy file
                        try:
                            with open(epotFile) as energyfile:
                                mosaics_energy = float(list(energyfile)[-1].split(' ')[0])
                            mosaics_energies[run_name] = [mosaics_energy, run_sequence, run_bppt]
                            mosaics_energy_cache[(run_sequence, run_bppt, remove_bases_A, remove_bases_B)] = mosaics_energy

                            os.remove(epotFile)
                            mosaics_energy = -1
                            
                        except (IndexError, FileNotFoundError) as _e:
                            mosaics_errorcount += 1
                            print('./pipeline.sh '+run_sequence+' '+str(run_bppt)+' '+str(remove_bases_A)+' '+str(remove_bases_B)+' '+inputFile+' '+pdbFile+' '+tempFile, "| MOSAICS not successful:", subprocess_output)
                            runs_success = False
                            break

                    else:
                        mosaics_energies[run_name] = None

                if runs_success:
                    currentOligo.add_energies(potIndelLength, potIndelStart, mosaics_energies)
                    progress_bar.update(1)
                
        data.append(currentOligo)

        try:
            os.remove(pdbFile)
            os.remove(inputFile)
            os.remove(epotFile)
        except FileNotFoundError as _e:
            pass


if __name__ == "__main__":
    import sys

    oligoFile  = "oligos.txt" if len(sys.argv) <= 1 else sys.argv[1]
    numWorkers = 48           if len(sys.argv) <= 2 else int(sys.argv[2])
    mosaicsWindow = 21        if len(sys.argv) <= 3 else int(sys.argv[3])
    indelFile  = ""           if len(sys.argv) <= 4 else sys.argv[4] # file containing true labels for indels, format like oligos.txt - not needed

    #TODO: test if this works on ARC when we don't give an oligo file
    # in caster.py: only calculate benchmarks for points that had a ground truth given
    # retrieve a requirements file from ARC environment
    
    startTime = time.time()

    global mosaics_energy_cache
    mosaics_energy_cache = None
    if os.path.exists("MOSAICS_EnergyCache.pickle"):
        energy_cache_file = open("MOSAICS_EnergyCache.pickle", "rb")
        mosaics_energy_cache = pickle.load(energy_cache_file)
        energy_cache_file.close()
        print("loaded energy cache")
    if not os.path.isdir("temp"): os.makedirs("temp")
    
    data, mosaics_errorcount = loadAllenData(indelFile, oligoFile, mosaicsWindow=mosaicsWindow, numWorkers=numWorkers)
    print("loading", len(data), "editing outcomes took", time.time()-startTime, "seconds, had to ignore", mosaics_errorcount, "sequences")
    
    pickle.dump(data, open("mosaics_data.pickle", "wb"))
    pickle.dump(mosaics_energy_cache, open("MOSAICS_EnergyCache.pickle", "wb"))
