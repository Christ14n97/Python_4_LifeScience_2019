# DNA Toolset/Code testing file

from bio_seq import bio_seq


DNAstr = "GCTGAAGCTTGGAACAGTGGAATCGCTATGCGTACGATGACGACGTGCAC"


test_dna = bio_seq()
test_dna.generate_rnd_seq(100, "RNA")

print(f'\n{test_dna.get_seq_info()}')

print(f'\nNucleotide Frequency: {test_dna.countNucFrequency()}')
print(f'\nGC%: {test_dna.gc_content()}%')
print(f'\nGC% in subseq: {test_dna.gc_content_subseq()}%')
print(f'\nAminoacids: {test_dna.translate_seq()}')
for rf in test_dna.gen_reading_frames():
    print(f'\n[ORF]: {rf}')

print(f'\nProteins: {test_dna.all_proteins_from_ORF()}')
