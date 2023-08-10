from random import sample
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from pathlib import Path


def random_permutation(seq):
    return ''.join(sample(list(seq), len(seq))).upper().replace('T', 'U')


def make_perm_record(og_header, perm_number, perm_seq):
    new_header = f'{og_header}_perm-{perm_number}'
    record = SeqRecord(Seq(perm_seq), id=new_header, description='')
    return record


def write_perm_record(perm_record, path):
    SeqIO.write(
        [perm_record], path, 'fasta'
    )


def main():

    seq = SeqIO.read(snakemake.input[0], 'fasta')
    bases = str(seq.seq).upper()
    num_perms = 100
    for i, each_path in enumerate(snakemake.output):
        perm = random_permutation(bases)
        perm_record = make_perm_record(seq.description, i, perm)
        write_perm_record(perm_record, each_path)


if __name__ == '__main__':
    main()
