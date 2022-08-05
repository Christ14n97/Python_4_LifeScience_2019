# DNA Toolkit file

# import this file 'from dnaToolkit import *'
# also 'import random'

from Structures import *
from collections import Counter
# check the sequence to make sure it is a DNA string


def validateSeq(dna_seq):
    """ check sequences to make sure it's a DNA string"""
    tmpseq = dna_seq.upper()
    for nuc in tmpseq:
        if nuc not in Nucleotides:
            return False
    return tmpseq

# generate a random str of DNA
    # import random
    # randDNAsrt=''.join([random.choice(Nucleotides) for nuc in range(50)])


def countNucFrequency(seq):
    """ return frequency of the nucleotides in a DNA string"""
    tmpFreqDict = {"A": 0, "C": 0, "G": 0, "T": 0}
    for nuc in seq:
        tmpFreqDict[nuc] += 1
    return tmpFreqDict


def transcription(seq):
    """ swapping timine with uracil and return a RNA string from a DNA string"""
    return seq.replace("T", "U")


def compleString(seq):
    """ return a complemetary string of a DNA single string"""
    return "".join([DNA_Complement[nuc] for nuc in seq])


def reverse_complement(seq):
    """
    Swapping adenine with thymine and guanine with cytosine.
    Reversing newly generated string
    """
    return ''.join([DNA_Complement[nuc] for nuc in seq])[::-1]

    # Pythonic approach. A little bit faster solution
    # mapping = str.maketrans('ATCG', 'TAGC')
    # return seq.translate(mapping)[::-1]


def gc_content(seq):
    """ GC content in a DNA/RNA seq. Usefull if we wanna see
    if a helix is thightly close (three bonds) or there's a interesting pattern"""

    return round(seq.count('C') + seq.count('G')/len(seq) * 100)


def gc_content_subseq(seq, k=5):
    """ GC content in DNA/RNA sub-sequence length k"""
    res = []
    for i in range(0, len(seq) - k + 1, k):
        subseq = seq[i:i + k]
        res.append(gc_content(subseq))
    return res


def translate_seq(seq, init_pos):
    """Translates a DNA sequence into an aa sequence"""
    return [DNA_Codons[seq[pos:pos+3]] for pos in range(init_pos, len(seq)-2, 3)]


def codon_usage(seq, aminoacid):
    """Provides the frequency of each codon encoding a given aminoacid in a DNA sequence"""
    tmpList = []
    for i in range(0, len(seq) - 2, 3):
        if DNA_Codons[seq[i:i + 3]] == aminoacid:
            tmpList.append(seq[i:i + 3])

    freqDict = dict(Counter(tmpList))
    # 'dict(Counter)'creates a dictionary with key:value where value is the frequence our key appears
    totalWight = sum(freqDict.values())
    for seq in freqDict:
        freqDict[seq] = round(freqDict[seq] / totalWight, 2)
    return freqDict


def gen_reading_frames(seq):
    """ Generate the six reading frames of a DNA sequence in aa, including reverse complement"""
    frames = []
    frames.append(translate_seq(seq, 0))
    frames.append(translate_seq(seq, 1))
    frames.append(translate_seq(seq, 2))
    frames.append(translate_seq(reverse_complement(seq), 0))
    frames.append(translate_seq(reverse_complement(seq), 1))
    frames.append(translate_seq(reverse_complement(seq), 2))
    return frames


def proteins_from_orf(aa_seq):
    """ Compute all possible proteins in an aminoacid seq and return a list of possible proteins"""
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


def all_proteins_from_ORF(seq, startReadPos=0, endReadPos=0, ordered=False):
    """ Compute all possible proteins for all ORF, which could be a specific range or from the entire seq"""
    """ Proteine search DB: https://www.ncbi.nlm.nih.gov/nuccore/NM_001185097.2"""
    """ API can be used to pull protein info"""
    if endReadPos > startReadPos:
        orf = gen_reading_frames(seq[startReadPos:endReadPos])
        # specifying the range
    else:
        orf = gen_reading_frames(seq)
        # from the entire seq
    res = []
    for rf in orf:
        prots = proteins_from_orf(rf)
        for p in prots:
            res.append(p)
    if ordered:
        return sorted(res, key=len, reverse=True)
    return res
