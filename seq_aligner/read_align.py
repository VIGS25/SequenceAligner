import json
from seq_aligner.genome_index import *
import numpy as np
from functools import partial
from operator import is_not


class readAlign(object):
    """ A class readAlign that creates the genome index, queries the genome, aligns
    the reads and produces the FLAG attribute for aligned reads"""

    def __init__(self, filename, params, method="kmer"):
        self.index_method = method
        self.filename = filename
        indexers = {'kmer':kmerIndex, 'suffixArray':suffixArray, 'FM':fmIndex}
        indexer = indexers[method]

        if method == "kmer":
            k = params.get('k', 16)
            kmer_spacing = params.get("kmer_spacing", 4)
            self.params = {"k":k, "kmer_spacing":kmer_spacing}
            self.indexer = indexer(k=k, kmer_spacing=kmer_spacing, genome_file=filename)

        elif method == "suffix_array":
            self.params = dict()
            self.indexer = indexer()

        else:
            raise NotImplementedError()

        self.indexer.generate_index()

    def get_read_index(self, read):
        """Computes the indices where a kmer of the read could be found."""
        # read: Read Sequence
        index = []

        if self.index_method == 'kmer':
            k = self.params['k']
            kmer_spacing = self.params['kmer_spacing']
            trials = range(kmer_spacing*5)

            forward_file, rev_comp_file = self.indexer.index_files

            with open(forward_file) as f:
                forward_index = json.load(f)
            with open(reverse_file) as f:
                rev_comp_index = json.load(f)

            for start_loc in trials:
                kmer = read[start_loc:start_loc+k]
                if kmer in forward_index or kmer in rev_index:
                    start = start_loc
                    break

            for i in range(start, len(read)-k + 1, kmer_spacing):
                index.append(read[i:i+k])

        elif self.index_mthod == 'FM':
            for i in range(0, len(read), len(read)//8):
                index.append(read[i: i + len(read)//8])

        return index

    def query_index(self, read):
        # read: Read Sequence
        # hitsFor: List of matches on forward strand
        # hitsRev: List of matches on reverse strand
        # index: Index of the forward strand
        # revIndex: Index of the reverse strand

        forward_hits = []
        reverse_hits = []

        if self.index_method == 'kmer':
            read_index = self.get_read_index(read)

            forward_file, rev_comp_file = self.indexer.index_files

            with open(forward_file) as f:
                forward_index = json.load(f)
            with open(rev_comp_file) as f:
                rev_comp_index = json.load(f)

            for kmer in read_index:
                if kmer in index:
                    hitsFor.append(index[kmer])
                elif kmer in revIndex:
                    hitsRev.append(revIndex[kmer])

        return hitsFor, hitsRev

    def getCIGAR(self, alignment, type='semiglobal'):
        """Computes the CIGAR string for the given alignment."""
        keys = ['M', 'I', 'D', 'N', 'S', 'H', 'P', '=', 'X']
        cigarCounts = {key: 0 for key in keys}

        # alignment is a list that contains both the read and genome alignment.
        cigar = ""

        if type == "semiglobal":
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

        else:
            raise NotImplementedError()


    def generate_edit_matrix(self, read, focus_genome, cost_method):

        cost_match = cost_method[0]
        cost_gap = cost_method[1]
        cost_mismatch = cost_method[2]

        edit_matrix = np.zeros((len(read) + 1, len(focus_genome) + 1), dtype='int64')
        edit_matrix[1:len(read) + 1, 0] = cost_gap * (np.arange(1, len(read) + 1))
        path_matrix = np.zeros(shape=edit_matrix.shape)

        for i in range(1, edit_matrix.shape[0]):
            for j in range(1, edit_matrix.shape[1]):

                ifEquals = [cost_match if read[i - 1] == focus_genome[j - 1] else cost_mismatch]
                edit_matrix[i, j] = min([edit_matrix[i - 1, j - 1] + ifEquals, edit_matrix[i - 1, j] + cost_gap,
                                        edit_matrix[i, j - 1] + cost_gap])
                pos = [edit_matrix[i - 1, j - 1] + ifEquals, edit_matrix[i - 1, j] + cost_gap,
                       edit_matrix[i, j - 1] + cost_gap].index(edit_matrix[i, j])
                if pos == 0:
                    path_matrix[i, j] = 1
                elif pos == 1:
                    path_matrix[i, j] = 2
                else:
                    path_matrix[i, j] = 3

        alignedRead = ''
        alignedGenome = ''
        minScore = np.min(edit_matrix[-1, :])

        possible_starts = int(np.argmin(edit_matrix[-1, :]))
        start = (len(read), possible_starts)

        next = start

        while next[0] != 0 and next[1] > 0:

            path = path_matrix[next[0], next[1]]

            if path == 1:

                alignedRead = read[next[0] - 1] + alignedRead
                alignedGenome = focus_genome[next[1]-1] + alignedGenome
                next = (next[0] - 1, next[1] - 1)

            elif path == 2:

                alignedRead = read[next[0] - 1] + alignedRead
                alignedGenome = '-' + alignedGenome
                next = (next[0] - 1, next[1])

            else:
                alignedRead = '-' + alignedRead
                alignedGenome = focus_genome[next[1]-1] + alignedGenome
                next = (next[0], next[1] - 1)

        return int(minScore), alignedRead, alignedGenome, next[1]

    def get_alignment_candidates(self, read, hits, repeat):

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


    def align(self, read, cost_method, method='semiglobal'):
        alignment = ['', '']
        pos = 0
        cigar = ''
        hits_forward, hits_reverse = self.query_index(read)

        if len(hitsF):
            cand_forward= self.get_alignment_candidates(read, hits_forward, repeat=len(hits_forward))

        if len(hitsR):
            cand_rev = self.get_alignment_candidates(read, hits_forward, repeat=len(hits_forward))

        possible_alignemnts = []

        if self.index_method == 'kmer':

            for hit in candidatesF:

                if hit > 1*len(read):

                    focus_genome = self.genome[hit - (1*len(read)): hit + (1*len(read))]
                    pos_start = hit - (1*len(read))

                else:
                    focus_genome = self.genome[0: hit + (1*len(read))]
                    pos_start = 0

                minScore, alignedRead, alignedGenome, startAlignment = self.generateedit_matrix(read, focus_genome, cost_method)

                possible_alignemnts.append((minScore, [alignedRead, 'F'+alignedGenome], pos_start + startAlignment + 1))

            for hit in candidatesR:

                revGenome = gI.reverseComplement(self.genome)

                if hit > 1*len(read):
                    focus_genome = revGenome[hit - (1*len(read)): hit+ (1*len(read))]
                    pos_start = hit - (1*len(read))
                else:
                    focus_genome = revGenome[0: hit + (2*len(read))]
                    pos_start = 0
                minScore, alignedRead, alignedGenome, startAlignment = self.generateedit_matrix(read, focus_genome, cost_method)

                possible_alignemnts.append((minScore, [alignedRead, 'R'+alignedGenome], pos_start + startAlignment + len(read)))

        bestAlignment = min(possible_alignemnts, key=lambda t: t[0])
        alignment = bestAlignment[1]

        if alignment[1][0] == 'R':
            pos = len(self.genome) - bestAlignment[2] + 1
        else:
            pos = bestAlignment[2]

        cigar = self.getCIGAR(alignment)

        return alignment, pos, cigar

    def get_FLAG_bits(self, strands_read_pair, first, num_align_read1=1, num_align_read2=1, paired_seq=True):
        """getFLAGbits gives us the bit values of the FLAG attribute of the read alignment output"""

        # Ignoring the last three bit flags mentioned in the documentation

        # strands_read_pair: Strand to which the reads from read_pair align, ['R', 'F'] or ['F', 'R']
        # num_align_read1: Number of alignments of read1 to the reference genome
        # num_align_read2: Number of alignments of read2 to the reference genome
        # first: If read1 is the first segment in template, default = True
        # paired_seq: If the sequencing experiment was paired end sequencing, default = True
        # possible_bits: Dictionary of bits as mentioned in the documentation

        possible_bits = {hex(1 << exponent): False for exponent in range(11)}

        if paired_seq:
            possible_bits['0x1'] = True

        # Check for mapping and update bits
        if num_align_read1 > 0 and num_align_read2 > 0:
            possible_bits['0x2'] = True

        if num_align_read2 == 0:
            possible_bits['0x4'] = True
            # Place a check for the assumptions listed in SAM_v1 file then

        if num_align_read1 > 0 and num_align_read2 == 0:
            possible_bits['0x8'] = True

        # Bits for reverse complements of current and next segment of template

        if strands_read_pair[0] == 'F' and strands_read_pair[1] == 'R':
            possible_bits['0x20'] = True
        elif strands_read_pair[0] == 'R' and strands_read_pair[1] == 'F':
            possible_bits['0x10'] = True

        # Bits for first and last segments of template

        if first:
            possible_bits['0x40'] = True
        else:
            possible_bits['0x80'] = True

        if num_align_read1 > 1:
            possible_bits['0x100'] = True

        return possible_bits
