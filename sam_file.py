import re


class Header(object):
    """Header class for SAM file. Contains attributes as defined in \
    the SAM file header specification."""

    def __init__(self, seq_name, genome_length, version=1.4, sorting_order="unsorted"):
        self._VN = version  # Format version
        self._SO = sorting_order # Sorting order
        self._SN = seq_name # Sequence name
        self._LN = genome_length # genome length

    @property
    def VN(self):
        return self._VN

    @property
    def SO(self):
        return self._SO

    @property
    def SN(self):
        return self._SN

    @property
    def LN(self):
        return self._LN

    def write(self, filename):
        with open(filename, 'w') as f:
            f.write('@HD' + '\t')
            f.write('VN:' + str(self.VN) + '\t')
            f.write('SO:' + str(self.SO) + '\n')

        with open(filename, 'a') as f:
            f.write('@SQ' + '\t')
            f.write('SN:' + str(self.SN) + '\t')
            f.write('LN:' + str(self.LN) + '\n')


class SAMLine(object):
    """A single record of the SAM File with all 11 mandatory fields."""

    def __init__(self):
        self._QNAME = '*' 
        self._FLAG = 0
        self._RNAME = '*'
        self._POS = 0
        self._MAPQ = 0
        self._CIGAR = '*'
        self._RNEXT = '*'
        self._PNEXT = 0
        self._TLEN = 0
        self._SEQ = '*'
        self._QUAL = '*'

    @property
    def QNAME(self):
        return self._QNAME

    @property
    def FLAG(self):
        return self._FLAG

    @property
    def RNAME(self):
        return self._RNAME

    @property
    def POS(self):
        return self._POS

    @property
    def MAPQ(self):
        return self._MAPQ

    @property
    def CIGAR(self):
        return self._CIGAR

    @property
    def RNEXT(self):
        return self._RNEXT

    @property
    def TLEN(self):
        return self._TLEN

    @property
    def SEQ(self):
        return self._SEQ

    @property
    def QUAL(self):
        return self._QUAL

    def set_base_fields(self, query_name, ref_name, 
        cigar, pos=0, map_quality=99, rnext='=', pnext=0, seq='', quality='*'):

        self._QNAME = query_name
        self._RNAME = ref_name
        self._POS = pos
        self._MAPQ = map_quality
        self._CIGAR = cigar
        self._RNEXT = rnext
        self._PNEXT = pnext
        self._SEQ = seq
        self._QUAL = quality

    def set_TLEN(self, pos, cigars):
        left_most_pos= min(pos)
        loc = pos.index(left_most_pos)
        right_most_pos = pos[len(pos)-loc-1] + sum(
            list(map(int, re.split('[MDIHSPX=]', cigars[len(pos) - loc - 1][:-1]))))

        diff = right_most_pos - left_most_pos
        if loc == 0:
            self._TLEN = diff
        else:
            self._TLEN = -diff

    def set_FLAG(self, flagcodes):
        self._FLAG = sum([int(i, 16) if flagcodes[i] else 0 for i in flagcodes])

    def write(self, filename):
        """Write the attributes to file, in append mode."""
        with open(filename, 'a') as f:
            f.write(self.QNAME + '\t')
            f.write(str(self.FLAG) + '\t')
            f.write(self.RNAME + '\t')
            f.write(str(self.POS) + '\t')
            f.write(str(self.MAPQ) + '\t')
            f.write(self.CIGAR + '\t')
            f.write(self.RNEXT + '\t')
            f.write(str(self.PNEXT) + '\t')
            f.write(str(self.TLEN) + '\t')
            f.write(self.SEQ + '\t')
            f.write(self.QUAL + '\n')