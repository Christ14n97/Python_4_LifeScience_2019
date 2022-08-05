from bio_structs import *
import random
from collections import Counter


class bio_seq:
    """DNA sequence class. Default value: ATCG, DNA, No label"""

    def __init__(self, seq="ATCG", seq_type="DNA", label="No Label"):
        """ Sequence initialization, validation"""
        self.seq = seq.upper()
        self.label = label
        self.seq_type = seq_type
        self.is_valid = self.__validateSeq()
        assert self.is_valid, f"Provided data does not seem to be correct {self.seq_type} sequence"
        # assert is the best way to test if something went correct, otherwise it will stop the execution

    def __validateSeq(self):
        """Check sequences to make sure if it is either a DNA or RNA string, retunrs true or false"""
        return set(Nucleotide_base[self.seq_type]).issubset(self.seq)

    def get_seq_biotype(self):
        """Returns sequence type"""
        return self.seq_type

    def get_seq_info(self):
        """Returns 4 strings. Full sequence information"""
        return f"[Label]: {self.label}\n[Sequence]: {self.seq}\n[Biotype]: {self.seq_type}\n[Length]: {len(self.seq)}"

    def generate_rnd_seq(self, length=10, seq_type="DNA"):
        """Generates a random DNA sequence, provided the length"""
        seq = ''.join([random.choice(Nucleotide_base[seq_type])
                       for nuc in range(length)])
        self.__init__(seq, seq_type, "Randomly generated sequence")
        # so, if we want to use a previous created variable, it's going to be checked

    def countNucFrequency(self):
        """Return frequency of the nucleotides in a DNA string, as a dictiorary"""
        return dict(Counter(self.seq))

    def transcription(self):
        """Swapping thymine with uracil and return a RNA string from a DNA string"""
        if self.seq_type == "DNA":
            return self.seq.replace("T", "U")
        return "Not a DNA sequence"

    def reverse_complement(self):
        """
        Swapping adenine with thymine and guanine with cytosine.
        Reversing newly generated string
        """
        if self.seq_type == "DNA":
            mapping = str.maketrans('ATCG', 'TAGC')
        else:
            mapping = str.maketrans('AUCG', 'UAGC')
        return self.seq.translate(mapping)[::-1]

    def gc_content(self):
        """GC content in a DNA/RNA seq. Usefull if we wannawant to check
        if a helix is thightly close (three bonds) or there's a interesting pattern"""
        return round(self.seq.count('C') + self.seq.count('G')/len(self.seq) * 100)

    def gc_content_subseq(self, k=5):
        """GC content in DNA/RNA sub-sequence length k"""
        res = []
        for i in range(0, len(self.seq) - k + 1, k):
            subseq = self.seq[i:i + k]
            res.append(round(subseq.count('C') +
                             subseq.count('G')/len(subseq) * 100))
        return res

    def translate_seq(self, init_pos=0):
        """Translates a DNA sequence into an amino acid sequence"""
        if self.seq_type == "DNA":
            return [DNA_Codons[self.seq[pos:pos+3]] for pos in range(init_pos, len(self.seq)-2, 3)]
        elif self.seq_type == "RNA":
            return [RNA_Codons[self.seq[pos:pos+3]] for pos in range(init_pos, len(self.seq)-2, 3)]

    def codon_usage(self, aminoacid):
        """Provides the frequency of each codon encoding a given aminoacid in a DNA sequence"""
        tmpList = []
        if self.seq_type == "DNA":
            for i in range(0, len(self.seq) - 2, 3):
                if DNA_Codons[self.seq[i:i + 3]] == aminoacid:
                    tmpList.append(self.seq[i:i + 3])
        elif self.seq_type == "RNA":
            for i in range(0, len(self.seq) - 2, 3):
                if RNA_Codons[self.seq[i:i + 3]] == aminoacid:
                    tmpList.append(self.seq[i:i + 3])

        freqDict = dict(Counter(tmpList))
        # 'dict(Counter)'creates a dictionary with key:value where value is the frequence our key appears
        totalWight = sum(freqDict.values())
        for self.seq in freqDict:
            freqDict[self.seq] = round(freqDict[self.seq] / totalWight, 2)
        return freqDict

    def gen_reading_frames(self):
        """Generate the six reading frames of a DNA sequence in aa, including reverse complement"""
        frames = []
        frames.append(self.translate_seq(0))
        frames.append(self.translate_seq(1))
        frames.append(self.translate_seq(2))
        tmp_seq = bio_seq(self.reverse_complement(), self.seq_type)
        frames.append(tmp_seq.translate_seq(0))
        frames.append(tmp_seq.translate_seq(1))
        frames.append(tmp_seq.translate_seq(2))
        del tmp_seq
        # so we can save some memory
        return frames

    def proteins_from_orf(self, aa_seq):
        """Compute all possible proteins in an aminoacid seq and return a list of possible proteins"""
        current_prot = []
        proteins = []

        for aa in aa_seq:
            if aa == "_":
                # STOP accumulating aa if '_' - STOP was found
                if current_prot:
                    for p in current_prot:
                        proteins.append(p)
                    current_prot = []
            else:
                # START accumulating aa if 'M' - START was found
                if aa == "M":
                    current_prot.append("")
                for i in range(len(current_prot)):
                    current_prot[i] += aa
        return proteins

    def all_proteins_from_ORF(self, startReadPos=0, endReadPos=0, ordered=False):
        """Compute all possible proteins for all ORF, which could be a specific range or from the entire seq"""
        """Proteine search DB: https://www.ncbi.nlm.nih.gov/nuccore/NM_001185097.2"""
        """API can be used to pull protein info"""

        if endReadPos > startReadPos:
            tmp_seq = bio_seq(self.seq[startReadPos:endReadPos], self.seq_type)
            orf = tmp_seq.gen_reading_frames()
            # specifying the range
        else:
            orf = self.gen_reading_frames()
            # from the entire seq

        res = []
        for rf in orf:
            prots = self.proteins_from_orf(rf)
            for p in prots:
                res.append(p)

        if ordered:
            return sorted(res, key=len, reverse=True)

        return res
