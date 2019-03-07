"""
Contains functions and classes to deal with genome indexing
"""

import os
import json

def get_complement(nucleotide):
    """Returns the complement of the given nucleotide."""
    complement_dict = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G', 'N': 'N'}
    return complement_dict[nucleotide]

def get_reverse_complement(sequence):
    """Returns the reverse complement of the given sequence."""
    complement = list(map(get_complement, sequence))
    rev_complement = list(reversed(complement))
    return "".join(rev_complement)

class Indexer(object):

    def __init__(self, **kwargs):
        self.genome, self.genome_name = self.read_genome(self.genome_file)
        self.index_files = [os.path.abspath("genomeIndex.txt"), os.path.abspath("revGenomeIndex.txt")]

    @staticmethod
    def read_genome(filename):
        """Reads the genome from the given filename."""
        genome = ""
        print("Reading the genome...")
        with open(filename, 'r') as f:
            genome_name = (f.readline()[1:].split(' ')[0].rstrip())
            for line in f:
                if not line[0] == '>':
                    genome += line.rstrip()

        print("Finished reading the genome." + "\n")
        return genome, genome_name

    @staticmethod
    def _build_fn(self, forward=True):
        raise NotImplementedError("Subclasses must implement for themselves.")

    def generate_index(self):
        print("[INFO]: Generating Index for the reference genome...")
        self._build_fn(forward=True)

        print("[INFO]: Generating Suffix Array Index for the reverse complement of \the reference genome...")
        self._build_fn(foward=False)


class kmerIndex(Indexer):
    """kmerIndex class for indexing the genome."""

    def __init__(self, k, kmer_spacing, genome_file, **kwargs):
        self.k = k
        self.kmer_spacing = kmer_spacing
        self.genome_file = genome_file
        super().__init__(**kwargs)

    def _build_fn(self, forward=True):
        index = dict()
        savefile = "forward_index.txt"
        genome = self.genome

        if not forward:
            savefile = "rev_complement_index.txt"
            genome = get_reverse_complement(self.genome)

        for i in range(0, len(genome) - self.k + 1, self.kmer_spacing):
            if genome[i:i + self.k] in index.keys():
                index[self.genome[i:i + self.k]].append(i)
            else:
                index[self.genome[i:i + self.k]] = [i]

        with open(os.path.abspath(savefile), "w") as f:
            f.write(json.dumps(index))


class suffixArray(Indexer):

    def __init__(self, **kwargs):
        super().__init__(**kwargs)

    def _build_fn(self, forward=True):
        """Builds the suffix array."""
        savefile = "forward_index.txt"
        genome = self.genome

        if not forward:
            savefile = "rev_complement_index.txt"
            genome = get_reverse_complement(self.genome)

        SATups = sorted([(genome[offset:], offset) for offset in range(len(genome))], key=lambda x: x[0])
        suffix_array = list(map(lambda x: x[1], SATups))

        with open(savefile, "w") as f:
            f.write(json.dumps(suffix_array))

'''FM-Index'''

def SA(s):
    sortedSuffixes = sorted([(s[i:], i) for i in range(len(s))], key=lambda x: x[0])
    sa = list(map(lambda x: x[1], sortedSuffixes))

    return sa


def bwtFromSA(t, sa=None):
    if sa is None:
        sa = SA(t)

    bwt = []

    for i in range(len(sa)):
        if sa[i] > 0:
            bwt.append(t[sa[i] - 1])
        else:
            bwt.append('$')

    return ''.join(i for i in bwt)


def fmCheckpoints(bwt, cpIval=4):

    checkPoints = {}
    tally = {}

    for c in bwt:
        if c not in tally:
            tally[c] = 0
            checkPoints[c] = []

    for i in range(0, len(bwt)):
        tally[bwt[i]] += 1
        if (i % cpIval) == 0:
            for c in tally.keys():
                checkPoints[c].append(tally[c])

    return checkPoints


class FmUtlities(object):

    @staticmethod
    def downsampleSuffixArray(sa, ssaIval=4):

        ssa = {}
        for i in range(len(sa)):

            if sa[i] % ssaIval == 0:
                ssa[i] = sa[i]
        return ssa

    def __init__(self, t, cpIval=4, ssaIval=4):

        if t[-1] != '$':
            t += '$'

        sa = SA(t)

        self.bwt = bwtFromSA(t, sa)
        self.indexFiles = []
        self.ssa = self.downsampleSuffixArray(sa, ssaIval)
        self.length = len(self.bwt)

        self.checkPoints = fmCheckpoints(self.bwt, cpIval)

        tots = {}
        for c in self.bwt:
            tots[c] = tots.get(c, 0) + 1

        self.first = {}
        totc = 0
        for c, count in sorted(tots.items()):
            self.first[c] = totc
            totc += count


class fmIndex(object):

    def __init__(self, cpIval=4, ssaIval=4):

        self.indexFiles = []
        self.cpIval = cpIval
        self.ssaIval = ssaIval
        self.genome = ''
        self.genomeName = ''

    def generateIndex(self, filename):

        self.genome, self.genomeName = readGenome(filename)

        fmUtils = FmUtlities(self.genome, self.cpIval, self.ssaIval)
        with open("BWT_genome.txt", "w") as f:
            f.write(fmUtils.bwt + '\n')
            f.write(json.dumps(fmUtils.first))

        self.indexFiles.append(os.path.abspath("BWT_genome.txt"))

        with open("ssa_genome.txt", "w") as f:
            f.write(json.dumps(fmUtils.ssa))

        self.indexFiles.append(os.path.abspath("ssa_genome.txt"))

        with open("cp_genome.txt", "w") as f:
            f.write(json.dumps(fmUtils.checkPoints))

        self.indexFiles.append(os.path.abspath("cp_genome.txt"))

        fmUtils = FmUtlities(reverseComplement(self.genome), self.cpIval, self.ssaIval)

        with open("BWT_revGenome.txt", "w") as f:
            f.write(fmUtils.bwt + '\n')
            f.write(json.dumps(fmUtils.first))

        self.indexFiles.append(os.path.abspath("BWT_revGenome.txt"))

        with open("ssa_revGenome.txt", "w") as f:
            f.write(json.dumps(fmUtils.ssa))

        self.indexFiles.append(os.path.abspath("ssa_revGenome.txt"))

        with open("cp_revGenome.txt", "w") as f:
            f.write(json.dumps(fmUtils.checkPoints))

        self.indexFiles.append(os.path.abspath("cp_revGenome.txt"))