
"""
bcr_models testing routines.

"""

import os

import bcr_models as bcr

# Example sequences
# -----------------

EX_HEAVY     = ('EVQLVESGPGLVQPGKSLRLSCVASGFTFSGYGMHWVRQAPGKGLEWIALIIYDESNKYYADSVKG'
                'RFTISRDNSKNTLYLQMSSLRAEDTAVFYCAKVKFYDPTAPNDYWGQGTLVTVSS')
EX_KAPPA     = ('DIVMTQTPSTLSASVGDRVTLTCKASQDISYLAWYQQKPGKAPKKLIYAASSLQSGVPSRFSGSGS'
                'GTDFTLTISSLQPEDFATYYCLQQNSNWTFGQGTKVDIK')
EX_ALPHA     = ('DSVIQMQGQVTFSENDSLFINCTYSTTGYPTLFWYVQYSGEGPQLLLQVTTANNKGSSRGFEATYD'
                'KGTTSFHLQKTSVQEIDSAVYYCAISDLSGGSNAKLAFGKGTKLSVK')
EX_BETA      = ('ITQTPKFLIGQEGQKLTLKCQQNFNHDTMYWYRQDSGKGLRLIYYSITENDLQKGDLSEGYDASRE'
                'KKSSFSLTVTSAQKNEMTVFLCASSIRLASAETLYFGSGTRLTVL')
NUCLEO_ALPHA = ('atggaaactctcctgggagtgtctttggtgattctatggcttcaactggctagggtgaac'
                'agtcaacagggagaagaggatcctcaggccttgagcatccaggagggtgaaaatgccacc'
                'atgaactgcagttacaaaactagtataaacaatttacagtggtatagacaaaattcaggt'
                'agaggccttgtccacctaattttaatacgttcaaatgaaagagagaaacacagtggaaga'
                'ttaagagtcacgcttgacacttccaagaaaagcagttccttgttgatcacggcttcccgg'
                'gcagcagacactgcttcttacttctgtgctacgggtcatagtggaggtagcaactataaa'
                'ctgacatttggaaaaggaactctcttaaccgtgaatccaaatatccagaaccctgac')
NUCLEO_BETA  = ('ATGGGCACCAGGCTCCTCTGCTGGGCAGCCCTGTGCCTCCTGGGGGCAGATCACACAGGT'
                'GCTGGAGTCTCCCAGACCCCCAGTAACAAGGTCACAGAGAAGGGAAAATATGTAGAGCTC'
                'AGGTGTGATCCAATTTCAGGTCATACTGCCCTTTACTGGTACCGACAAAGCCTGGGGCAG'
                'GGCCCAGAGTTTCTAATTTACTTCCAAGGCACGGGTGCGGCAGATGACTCAGGGCTGCCC'
                'AACGATCGGTTCTTTGCAGTCAGGCCTGAGGGATCCGTCTCTACTCTGAAGATCCAGCGC'
                'ACAGAGCGGGGGGACTCAGCCGTGTATCTCTGTGCCAGCAGCTTAAGCGGGAGGGCTTGG'
                'ACAGATACGCAGTATTTTGGCCCAGGCACCCGGCTGACAGTGCTCGAGGACCTGAAAAAC'
                'GTG')

# Files
# -----

TESTDIR   = os.path.dirname(__file__)
KAPPA_HMM = os.path.join(TESTDIR, 'kappa.hmm')
HEAVY_HMM = os.path.join(TESTDIR, 'heavy.hmm')


# Pseudo classes for testing
# --------------------------

def IgComplex():
    pass


def IgChain_kappa():
    ig = bcr.IgChain(EX_KAPPA, TestTemplateDatabase(), TestPDBDatabase())
    ig.hmmalign(KAPPA_HMM)
    ig.canonical_structures(TestCsDatabase())
    return ig


def IgChain_heavy():
    ig = bcr.IgChain(EX_HEAVY, TestTemplateDatabase(), TestPDBDatabase())
    ig.hmmalign(HEAVY_HMM)
    ig.canonical_structures(TestCsDatabase())
    return ig


def TestCsDatabase():
    filename = os.path.join(TESTDIR, 'canonical_structures.json')
    with open(filename) as f:
        return bcr.db.CsJSONDatabase(f)


def TestPDBDatabase():
    return bcr.db.PDBDirectoryDatabase(os.path.join(TESTDIR, 'pdb'))


def TestTemplateDatabase():
    with open(os.path.join(TESTDIR, 'templates.csv')) as f:
        return bcr.db.TemplateCSVDatabase(f)


def run_tests():
    import doctest
    import unittest

    import bcr_models.modelling
    import bcr_models.utils

    suite = unittest.TestSuite()
    suite.addTests(doctest.DocTestSuite(bcr.db))
    suite.addTests(doctest.DocTestSuite(bcr.utils))
    suite.addTests(doctest.DocTestSuite(bcr.modelling))
    suite.addTests(doctest.DocTestSuite(bcr.canonical_structures))

    runner = unittest.TextTestRunner()
    runner.run(suite )


if __name__ == '__main__':
    run_tests()