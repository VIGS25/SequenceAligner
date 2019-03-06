import json
import genomeIndexing as gI
import numpy as np
from functools import partial
from operator import is_not


class readAlign(object):
    """ A class readAlign that creates the genome index, queries the genome, aligns
    the reads and produces the FLAG attribute for aligned reads"""

    def __init__(self):
        # indexFiles: Contains the path to the index files
        # genome: Contains the reference genome
        # k: kmer size (useful only if kmer index is used)

        self.indexFiles = []
        self.genome = ''
        self.indexMethod = ''
        self.params = []
        self.genomeName = ''
        self.genomeLength = 0

    def createGenomeIndex(self, filename, params, method = 'kmer'):
        # Add an error to ensure correcting filename path
        self.indexMethod = method
        self.params = params

        if self.indexMethod == 'kmer':
            kmerIndex = gI.kmerIndex(k=self.params[0], kmerSpacing=self.params[1])
            kmerIndex.generateIndex(filename)
            self.genome, self.genomeName = kmerIndex.genome, kmerIndex.genomeName
            self.genomeLength = len(self.genome)
            self.indexFiles = kmerIndex.indexFiles

        elif self.indexMethod == 'suffixArray':
            SA = gI.suffixArray()
            SA.generateIndex(filename)
            self.genome, self.genomeName = SA.genome, SA.genomeName
            self.genomeLength = len(self.genome)
            self.indexFiles = SA.indexFiles

        elif self.indexMethod == 'FM':
            fmIndex = gI.fmIndex(cpIval=self.params[0], ssaIval=self.params[1])
            fmIndex.generateIndex(filename)
            self.genome, self.genomeName = fmIndex.genome, fmIndex.genomeName
            self.genomeLength = len(self.genome)
            self.indexFiles = fmIndex.indexFiles

    def getReadIndex(self, read):
        # read: Read Sequence
        # kmers: List of kmers in the read
        index = []

        if self.indexMethod == 'kmer':
            k = self.params[0]
            kmerSpacing = self.params[1]
            trials = range(kmerSpacing*5)

            with open(self.indexFiles[0]) as f:
                forIndex = json.load(f)
            with open(self.indexFiles[1]) as f:
                revIndex = json.load(f)

            for startLoc in trials:
                kmer = read[startLoc:startLoc+k]
                if kmer in forIndex or kmer in revIndex:
                    start = startLoc
                    break

            for i in range(start, len(read)-k + 1, kmerSpacing):
                index.append(read[i:i+k])

        elif self.indexMethod == 'FM':
            for i in range(0, len(read), len(read)//8):
                index.append(read[i: i + len(read)//8])

        return index

    def queryIndex(self, read):
        # read: Read Sequence
        # hitsFor: List of matches on forward strand
        # hitsRev: List of matches on reverse strand
        # index: Index of the forward strand
        # revIndex: Index of the reverse strand

        hitsFor = []
        hitsRev = []

        if self.indexMethod == 'kmer':
            readIndex = self.getReadIndex(read)

            with open(self.indexFiles[0]) as f:
                index = json.load(f)
            with open(self.indexFiles[1]) as f:
                revIndex = json.load(f)

            for kmer in readIndex:
                if kmer in index:
                    hitsFor.append(index[kmer])
                elif kmer in revIndex:
                    hitsRev.append(revIndex[kmer])

        elif self.indexMethod == 'FM':
            readIndex = self.getReadIndex(read)

            with open(self.indexFiles[0], 'r') as f:
                bwtF_L = f.readline().rstrip()
                bwtF_F = json.load(f)

            with open(self.indexFiles[1], 'r') as f:
                smallSA_F = json.load(f)

            with open(self.indexFiles[2], 'r') as f:
                checkPointsF = json.load(f)

            with open(self.indexFiles[3], 'r') as f:
                bwtR_L = f.readline().rstrip()
                bwtR_F = json.load(f)

            with open(self.indexFiles[4], 'r') as f:
                smallSA_R = json.load(f)

            with open(self.indexFiles[5], 'r') as f:
                checkPointsR = json.load(f)

            def rank(bwt, c, row, checkPoints, cpIval):
                if row < 0 or c not in checkPoints:
                    return 0
                i, nocc = row, 0
                # Always walk to left (up) when calculating rank
                while (i % cpIval) != 0:
                    if bwt[i] == c:
                        nocc += 1
                    i -= 1
                return checkPoints[c][i // cpIval] + nocc

            def count(c, BWT_F):
                if c not in BWT_F:
                    # (Unusual) case where c does not occur in text

                    for cc in sorted(BWT_F.keys()):
                        if c < cc:
                            return BWT_F[cc]

                    return BWT_F[cc]
                else:
                    return BWT_F[c]

            def getRange(p, BWT_L, BWT_F, checkPoints, cpIval):
                l, r = 0, len(BWT_L) - 1
                for i in range(len(p) - 1, -1, -1):
                    l = rank(BWT_L, p[i], l - 1, checkPoints, cpIval) + count(p[i], BWT_F)
                    r = rank(BWT_L, p[i], r, checkPoints, cpIval) + count(p[i], BWT_F) - 1
                    if r < l:
                        break
                return l, r + 1

            def resolve(row, BWT_L, BWT_F, checkPoints, cpIval, ssa):
                def stepLeft(row, BWT_L, BWT_F, checkPoints, cpIval):

                    c = BWT_L[row]
                    return rank(BWT_L, c, row - 1, checkPoints, cpIval) + count(c, BWT_F)

                nsteps = 0
                while row not in ssa:
                    row = stepLeft(row, BWT_L, BWT_F, checkPoints, cpIval)
                    nsteps += 1
                return ssa[row] + nsteps

            def hasSubstring(p, BWT_L, BWT_F, checkPoints, cpIval):
                l, r = getRange(p, BWT_L, BWT_F, checkPoints, cpIval)
                return r > l

            def occurrences(p, BWT_L, BWT_F, checkPoints, cpIval, ssa):
                l, r = getRange(p, BWT_L, BWT_F, checkPoints, cpIval)
                return [resolve(x, BWT_L, BWT_F, checkPoints, cpIval, ssa) for x in range(l, r)]

        return hitsFor, hitsRev

    def getCIGAR(self, alignment, type='semiglobal'):

        keys = ['M', 'I', 'D', 'N', 'S', 'H', 'P', '=', 'X']
        cigarCounts = {key: 0 for key in keys}

        # alignments is a list that contains both the read and genome alignment.
        cigar = ''

        if type == 'semiglobal':

            pos = 0

            while pos < len(alignment[0]):

                while alignment[0][pos] != ' ' and alignment[1][pos] != ' ':
                    cigarCounts['M'] += 1
                    pos += 1
                    if pos == len(alignment[0]):
                        break

                cigar += (str(cigarCounts['M']) + 'M')
                cigarCounts['M'] = 0

                if pos != len(alignment[0]):

                    if alignment[1][pos] == ' ':
                        cigar += str(cigarCounts['I'] + 1) + 'I'
                        pos += 1
                    else:
                        cigar += str(cigarCounts['D'] + 1) + 'D'
                        pos += 1

        return cigar

    def generateEditMatrix(self, read, focusGenome, costMethod):

        cost_match = costMethod[0]
        cost_gap = costMethod[1]
        cost_mismatch = costMethod[2]

        editMatrix = np.zeros((len(read) + 1, len(focusGenome) + 1), dtype='int64')
        editMatrix[1:len(read) + 1, 0] = cost_gap * (np.arange(1, len(read) + 1))
        pathMatrix = np.zeros(shape=editMatrix.shape)

        for i in range(1, editMatrix.shape[0]):
            for j in range(1, editMatrix.shape[1]):

                ifEquals = [cost_match if read[i - 1] == focusGenome[j - 1] else cost_mismatch]
                editMatrix[i, j] = min([editMatrix[i - 1, j - 1] + ifEquals, editMatrix[i - 1, j] + cost_gap,
                                        editMatrix[i, j - 1] + cost_gap])
                pos = [editMatrix[i - 1, j - 1] + ifEquals, editMatrix[i - 1, j] + cost_gap,
                       editMatrix[i, j - 1] + cost_gap].index(editMatrix[i, j])
                if pos == 0:
                    pathMatrix[i, j] = 1
                elif pos == 1:
                    pathMatrix[i, j] = 2
                else:
                    pathMatrix[i, j] = 3

        alignedRead = ''
        alignedGenome = ''
        minScore = np.min(editMatrix[-1, :])

        possible_starts = int(np.argmin(editMatrix[-1, :]))
        start = (len(read), possible_starts)

        next = start

        while next[0] != 0 and next[1] > 0:

            path = pathMatrix[next[0], next[1]]

            if path == 1:

                alignedRead = read[next[0] - 1] + alignedRead
                alignedGenome = focusGenome[next[1]-1] + alignedGenome
                next = (next[0] - 1, next[1] - 1)

            elif path == 2:

                alignedRead = read[next[0] - 1] + alignedRead
                alignedGenome = '-' + alignedGenome
                next = (next[0] - 1, next[1])

            else:
                alignedRead = '-' + alignedRead
                alignedGenome = focusGenome[next[1]-1] + alignedGenome
                next = (next[0], next[1] - 1)

        return int(minScore), alignedRead, alignedGenome, next[1]

    def getAlignmentCandidates(self, read, hits, repeat):

        candidates = []
        limit = 0

        hits = np.array(hits)
        candidates.extend(hits[0])

        while limit < min(repeat, 4):
            focus = hits[limit]

            for i in range(len(focus)):
                for j in range(limit + 1, len(hits)):

                    pos = list(map(lambda x: x if abs(x - focus[i]) > len(read) else None, hits[j]))
                    candidates.extend(pos)

            limit += 1

        candidates = list(filter(partial(is_not, None), candidates))
        return list(set(candidates))


    def align(self, read, costMethod, method='semiglobal'):

        alignment = ['', '']
        pos = 0
        cigar = ''

        hitsF, hitsR = self.queryIndex(read)
        candidatesF = []
        candidatesR = []

        if len(hitsF):
            candidatesF = self.getAlignmentCandidates(read, hitsF, repeat=len(hitsF))

        if len(hitsR):
            candidatesR = self.getAlignmentCandidates(read, hitsR, repeat=len(hitsR))

        possibleAlignments = []

        if self.indexMethod == 'kmer':

            for hit in candidatesF:

                if hit > 1*len(read):

                    focusGenome = self.genome[hit - (1*len(read)): hit + (1*len(read))]
                    pos_start = hit - (1*len(read))

                else:
                    focusGenome = self.genome[0: hit + (1*len(read))]
                    pos_start = 0

                minScore, alignedRead, alignedGenome, startAlignment = self.generateEditMatrix(read, focusGenome, costMethod)

                possibleAlignments.append((minScore, [alignedRead, 'F'+alignedGenome], pos_start + startAlignment + 1))

            for hit in candidatesR:

                revGenome = gI.reverseComplement(self.genome)

                if hit > 1*len(read):
                    focusGenome = revGenome[hit - (1*len(read)): hit+ (1*len(read))]
                    pos_start = hit - (1*len(read))
                else:
                    focusGenome = revGenome[0: hit + (2*len(read))]
                    pos_start = 0
                minScore, alignedRead, alignedGenome, startAlignment = self.generateEditMatrix(read, focusGenome, costMethod)

                possibleAlignments.append((minScore, [alignedRead, 'R'+alignedGenome], pos_start + startAlignment + len(read)))

        bestAlignment = min(possibleAlignments, key=lambda t: t[0])
        alignment = bestAlignment[1]

        if alignment[1][0] == 'R':
            pos = len(self.genome) - bestAlignment[2] + 1
        else:
            pos = bestAlignment[2]

        cigar = self.getCIGAR(alignment)

        return alignment, pos, cigar

    def getFLAGbits(self, strandsReadPair, first, nAlignRead1=1, nAlignRead2=1, pairedSeq = True):

        """getFLAGbits gives us the bit values of the FLAG attribute of the read alignment output"""

        # Ignoring the last three bit flags mentioned in the documentation

        # strandsReadPair: Strand to which the reads from readPair align, ['R', 'F'] or ['F', 'R']
        # nAlignRead1: Number of alignments of read1 to the reference genome
        # nAlignRead2: Number of alignments of read2 to the reference genome
        # first: If read1 is the segment in template, default = True
        # pairedSeq: If the sequencing experiment was paired end sequencing, default = True
        # possible_bits: Dictionary of bits as mentioned in the documentation

        possible_bits = {hex(1 << exponent): False for exponent in range(11)}

        if pairedSeq:
            possible_bits['0x1'] = True

        # Check for mapping and update bits
        if nAlignRead1 > 0 and nAlignRead2 > 0:
            possible_bits['0x2'] = True

        if nAlignRead1 == 0:
            possible_bits['0x4'] = True
            # Place a check for the assumptions listed in SAM_v1 file then

        if nAlignRead1 > 0 and nAlignRead2 == 0:
            possible_bits['0x8'] = True

        # Bits for reverse complements of current and next segment of template

        if strandsReadPair[0] == 'F' and strandsReadPair[1] == 'R':
            possible_bits['0x20'] = True
        elif strandsReadPair[0] == 'R' and strandsReadPair[1] == 'F':
            possible_bits['0x10'] = True

        # Bits for first and last segments of template

        if first:
            possible_bits['0x40'] = True
        else:
            possible_bits['0x80'] = True

        if nAlignRead1 > 1:
            possible_bits['0x100'] = True

        return possible_bits











