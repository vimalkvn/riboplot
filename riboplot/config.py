# -*- coding: utf-8 -*-
import os


class Config(object):
    # get the running directory of this file, move one level up to get the
    # application directory
    APP_DIR = os.path.split(os.path.dirname(os.path.abspath(__file__)))[0]
    PKG_DATA_DIR = os.path.join(APP_DIR, 'riboplot', 'data')


class TestingConfig(Config):
    """Testing configuration"""
    # directory for test data sets
    TEST_DATA_DIR = os.path.join(Config.APP_DIR, 'tests/data')

    # sorted BAM format Ribo-Seq alignment file
    RIBO_FILE = os.path.join(TEST_DATA_DIR, '5hRPFsorted.bam')

    # sorted BAM format RNA-Seq alignment file
    RNA_FILE = os.path.join(TEST_DATA_DIR, '5hmRNAsorted.bam')

    # name of the transcript to use (riboplot)
    TRANSCRIPT_NAME = 'gi|148357119|ref|NM_001098396.1|'

    # length of the above transcript (calculated manually)
    TRANSCRIPT_LENGTH = 1298

    # transcript sequence
    TRANSCRIPT_SEQUENCE = (
        'ATAGGACTCCGGTTCTATTTTGTGGGTTTCTGGAACCCGGGGCCATGATTGAGAGGGACGGCCGGGGGCATTCGTATTGC'
        'GCCGCTAGAGGTGAAATTCTTGGACCGGCGCAAGACGGACGAAAGCGAAAGCATTTGCCAAGAATGTTTTCATTAATCAA'
        'GAACGAAAGTCGGAGGTTCGAAGACGATCAGATACCGTCGTAGTTCCGACCGTAAACGATGCCGACCCGCGATCCGGCGG'
        'CGTTATTCCCATGACCCGCCGGGCAGCGTGCGGGAAACCACGAGTCTTTGGGTTCCGGGGGGAGTATGGTTGCAAAGCTG'
        'AAACTTAAAGGAATTGACGGAAGGGCACCACCAGGAGTGGAGCCTGCGGCTTAATTTGACTCAACACGGGAAACCTCACC'
        'CGGCCCGGACACGGAAAGGATTGACAGATTGATAGCTCTTTCTCGATTCTGTGGGTGGTGGTGCATGGCCGTTCTTAGTT'
        'GGTGGAGCGATTTGTCTGGTTCATTCCGATAACGAACGAGACTCCGGCATGCTAAATAGTTACGCGGCCCCGCGCGGTCG'
        'GCGTCCCAACTTCTTAGAGGGACAAGTGGCGTTCAGCCACGCGAGATGGAGCAATAACAGGTCTGTGATGCCCTTAGATG'
        'TCCGGGGCTGCACGCGCGCCACAATGGGCGGATCAACGTGTGCCTACCCTGCGCCGAGAGGCGCGGGTAACCCGTTGAAC'
        'CCCGCTCGTGATTGGGACTGGGGCTTGAAACTGTTTCCCATCAACGAGGAATTCCCAGTAAGCGCAGGTCATAAGCTTGC'
        'GTTGATTAAGTCCCTGCCCTTTGTACACACCGCCCGTCGCTACTACCGATTGGATGGTTTAGTGAGGTCCTCGGATCGGC'
        'CCCGCCGGGGCTCCTCGCCGGGCCCTGGCGGAGCGCCGAGAAGACGATCAAACTTGATCCTCTAGAGGAAGTAAAAGTCG'
        'TAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTACGGGTTCGAGAGGGGATCTCCCCCCCGCTCATTCCTCCGC'
        'TTCGGGGTCCCTCCGGGGTTACCCCAGGCTCGGAAAACGGTGAACCTGGCGCGGTCGGCGACAGCGAGCTCCCCGACGGG'
        'TACCCGCCTGGCCTCTCGGGGCCGTGGGTTCAAAGACCTTCCCGTCTCGACGGGGAGGCGCCCGTCCGGGGCTCTGCCGA'
        'CCGCCCGTTTTTAATAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA'
        'AAAAAAAAAAAAAAAAAA')

    # fasta file of the transcriptome
    TRANSCRIPTOME_FASTA = os.path.join(TEST_DATA_DIR, 'zebrafish.fna')

    # unrelated fasta file
    UNRELATED_FASTA = os.path.join(TEST_DATA_DIR, 'unrelated.fna')

    # unrelated BAM file
    UNRELATED_BAM = os.path.join(TEST_DATA_DIR, 'unrelated.bam')


class ProductionConfig(Config):
    """Production configuration"""
    # Additional variables can be listed here
    pass
