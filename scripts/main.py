#!/usr/bin/env python
import argparse
import re
import subprocess
import abc
from Bio import SeqIO

__author__ = 'Sergey Knyazev'


class Chain:
    __metaclass__ = abc.ABCMeta

    def __init__(self, query):
        self.__query = ''
        self._cdr1_begin = 0
        self._cdr1_length = 0
        self._cdr2_begin = 0
        self._cdr2_length = 0
        self._cdr3_begin = 0
        self._cdr3_length = 0
        self.query = query

    @abc.abstractclassmethod
    def _annotate(self):
        raise NotImplementedError

    @property
    def query(self):
        return self.__query

    @query.setter
    def query(self, query):
        self.__query = query
        self._annotate()

    def fr1(self):
        return self.query[:self._cdr1_begin-1]

    def fr2(self):
        return self.query[self._cdr1_begin+self._cdr1_length+1:self._cdr1_begin-1]

    def fr3(self):
        return self.query[self._cdr2_begin+self._cdr2_length+1:]

    def cdr1(self):
        return self.query[self._cdr1_begin:self._cdr1_begin+self._cdr1_length]

    def cdr1_begin(self):
        return self._cdr1_begin

    def cdr1_length(self):
        return self._cdr1_length

    def cdr2(self):
        return self.query[self._cdr2_begin:self._cdr2_begin+self._cdr2_length]

    def cdr2_begin(self):
        return self._cdr2_begin

    def cdr2_length(self):
        return self._cdr2_length

    def cdr3(self):
        return self.query[self._cdr3_begin:self._cdr3_begin+self._cdr3_length]

    def cdr3_begin(self):
        return self._cdr3_begin

    def cdr3_length(self):
        return self._cdr3_length


class HeavyChain(Chain):
    def __init__(self, query):
        super(HeavyChain, self).__init__(query)

    def _annotate(self):
        #find CDR1
        p = re.compile('C[A-Z]{13,15}W[IVFYAMLN][RKQVNC][QKHELR]')
        r = p.search(self.query)

        if r is None:
            print("CDR1 is not found in: " + self.query)
            exit(1)

        self._cdr1_begin = r.start() + 9
        self._cdr1_length = r.end() - self._cdr1_begin - 4

        #find CDR3
        p = re.compile('C[A-Z]{5,27}WG[A-Z][GRDS]')
        r = p.search(self.query, self._cdr1_begin + self._cdr1_length)

        if r is None:
            print("CDR3 is not found in: " + self.query)
            exit(1)

        self._cdr3_begin = r.start() + 3
        self._cdr3_length = r.end() - self._cdr3_begin - 4

        self._cdr2_begin = self._cdr1_begin + self._cdr1_length + 15
        self._cdr2_length = self._cdr3_begin - self._cdr2_begin - 33

        if self._cdr2_length < 4:
            print("CDR2 is not found in: " + self.query)
            exit(1)


#TODO: Implement class LightChain
class LightChain(Chain):
    def _annotate(self):
        pass


# TODO:
def find_homolog(sequence, blast_parameters):
    blast_results = subprocess.Popen(['blastlp', '-p', 'blastp',
                                      '-d', blast_parameters.database_path + '/' + blast_parameters.database_name,
                                      '-e', str(blast_parameters.e_value), '-M', blast_parameters.matrix,
                                      '-W', str(blast_parameters.word_size), '-m', str(blast_parameters.output_type),
                                      '-b', str(blast_parameters.number_to_keep)], stdout=subprocess.PIPE,
                                     stdin=subprocess.PIPE, stderr=subprocess.STDOUT).communicate(input=sequence)
    open('test.txt', mode='w').write(blast_results)


class BlastParameters():
    def __init__(self, database_path):
        self.database_name = ''
        self.database_path = database_path
        self.e_value = 0.0
        self.matrix = ''
        self.word_size = -1
        self.number_to_keep = 600
        self.output_type = 9


# TODO:
def find_fragments_homologs(chain, args):
    blast_parameters = BlastParameters(args.database)
    blast_parameters.e_value = 0.00001
    blast_parameters.matrix = "BLOSUM62"
    blast_parameters.word_size = 0
    blast_parameters.database_name = "hfr"
    fr_homolog = find_homolog(chain.fr1() + chain.fr2() + chain.fr3(), blast_parameters)
    blast_parameters.e_value = 2000
    blast_parameters.matrix = "PAM30"
    blast_parameters.word_size = 2
    blast_parameters.database_name = "h1"
    cdr1_homolog = find_homolog(chain.cdr1(), blast_parameters)
    blast_parameters.database_name = "h2"
    cdr2_homolog = find_homolog(chain.cdr2(), blast_parameters)
    blast_parameters.database_name = "h3"
    cdr3_homolog = find_homolog(chain.cdr3(), blast_parameters)
    return ''


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--input-fasta', help='input fasta file', required=True)
    parser.add_argument('-o', '--out-pdb', help='out pdb file', required=True)
    parser.add_argument('-d', '--database', help='database folder', required=True)
    parser.add_argument('-c', '--camelid', help='if query is camelid')
    return parser.parse_args()


def main():
    args = parse_args()
    heavy_chain = HeavyChain(str(SeqIO.parse(args.input_fasta, "fasta").__next__().seq))
    fragments_homologs = find_fragments_homologs(heavy_chain, args)

if __name__ == "__main__":
    main()

