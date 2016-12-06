
"""
Canonical structures
"""

import os
import sys
import json
import pickle

import numpy as np

import bcr_models as igm


try:
    FileNotFoundError
except NameError:
    FileNotFoundError = IOError


class CsDatabase(object):
    """
    Database interface for canonical structures.
    """

    def get(self, ig_chain):
        return {cdr: 0 for cdr in ig_chain.cdr_regions}

    def __repr__(self):
        return 'CsDatabase()'


class CsJSONDatabase(CsDatabase):
    """
    Canonical structure database from JSON file.
    """

    def __init__(self, jsonhandle):
        try:
            self._filename = jsonhandle.name
        except AttributeError:
            self._filename = 'Unknown'

        self._parse_json(jsonhandle)

    def get(self, ig_chain):
        """Calculate canonical structure from a cdr.

        :return: CDR cs types.
        :rtype:  dict

        """
        cs = {}
        for cdr, (s, e) in ig_chain.cdr_regions.items():
            cs[cdr] = 0

            #Go through the different possible cs
            cdrlen = str(len(ig_chain.aligned_seq[s:e].replace('-', '')))
            opts = self._db.get(cdr, {}).get(cdrlen, {'': 0})

            for opt, p_cs in opts.items():
                for pos, res in [o.split(';') for o in opt.split()]:
                    i = ig_chain.renum(pos)
                    if (ig_chain.aligned_seq[i] not in res and
                        ig_chain.leftaln_seq[i] not in res):
                        break
                else:
                    cs[cdr] = p_cs
                    break

        return cs

    def _parse_json(self, jsonhandle):
        """Parse the json filehandle."""
        self._db = json.load(jsonhandle)

    def __str__(self):
        return 'CsJSONDatabase({})'.format(self._filename)


class BuiltinCsDatabase(CsJSONDatabase):
    """
    Built in canonical structures.
    """

    def __init__(self):
        lcl = os.path.dirname(__file__)
        self._filename = os.path.join(lcl, 'data', 'canonical_structures.json')

        with open(self._filename) as f:
            self._parse_json(f)

        self.rfdb = RandomForestCsDatabase()

    def get(self, ig_chain):
        """Return RF database prediction if available, fall back to JSON tree."""
        cs_forest = self.rfdb.get(ig_chain)
        cs_tree = super(BuiltinCsDatabase, self).get(ig_chain)
        for cdr, cs_type in cs_forest.items():
            if not cs_type:
                cs_forest[cdr] = cs_tree[cdr]

        return cs_forest


class RandomForestCsDatabase(CsDatabase):
    """
    Predict canonical_structures using random forests.

    """

    def __init__(self, rf_dir=None):
        lcl = os.path.dirname(__file__)
        self._rf_dir = rf_dir
        if not self._rf_dir:
            self._rf_dir = os.path.join(lcl, 'data', 'rf{}'.format(sys.version_info.major))
        self._rf = {}

    def get(self, ig_chain):
        cs = {}
        seq = igm.utils.encode_sparse(ig_chain.aligned_seq)
        for cdr in ig_chain.cdr_regions:
            ln = ig_chain.cdr_len(cdr)
            rf = self.get_rf(cdr)
            cs[cdr] = rf.predict(seq + [ln])[0]
            if max(rf.predict_proba(seq + [ln])[0]) < 0.25:
                cs[cdr] = 0

        return cs

    def get_rf(self, cdr):
        """Load an RF from file.

        Usage::

            >>> csdb = RandomForestCsDatabase()
            >>> csdb.get_rf('A1') #doctest: +ELLIPSIS
            RandomForestClassifier(...)

        Fallback::

            >>> csdb.get_rf('Z') #doctest: +ELLIPSIS
            DummyRF()

        """

        if cdr not in self._rf:
            try:
                with open(os.path.join(self._rf_dir, '{}.pickle'.format(cdr)), 'rb') as f:
                    self._rf[cdr] = pickle.load(f)
            except FileNotFoundError:
                self._rf[cdr] = DummyRF()

        return self._rf[cdr]


class DummyRF:
    """
    Dummy predictor always returning 0
    """

    def predict(self, X):
        return [0, ]

    def predict_proba(self, X):
        return [[0, ]]

    def __bool__(self):
        return False

    def __repr__(self):
        return 'DummyRF()'


def train_rf_database(entries, cdr, hmm=None, n_estimators=100):
    """Train a random forest databse for canonical_structures."""
    data = []
    clusters = []
    for entry in entries:
        seq = entry['sequence'].replace('-', '')
        ig_chain = igm.IgChain(seq)
        if not hmm:
            ig_chain.hmmsearch(*igm.db.builtin_hmms())
            hmm = ig_chain._hmm
        else:
            ig_chain.hmmalign(hmm)

        datum = igm.utils.encode_sparse(ig_chain.aligned_seq)
        datum.append(ig_chain.cdr_len(cdr))
        data.append(datum)
        clusters.append(entry[cdr])

    data = np.array(data)

    import sklearn.ensemble
    rf = sklearn.ensemble.RandomForestClassifier(n_estimators=n_estimators, oob_score=True)
    rf.fit(data, clusters)

    return rf

