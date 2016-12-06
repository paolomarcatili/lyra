#!/usr/bin/env python

try:
    import cStringIO as StringIO
except ImportError:
    import io as StringIO

from flask import Flask
from flask import jsonify
from flask import request
from flask import render_template

import bcr_models as bcr

app = Flask(__name__)

hmms = bcr.db.builtin_hmms()
template_db = bcr.db.BuiltinTemplateDatabase()
pdb_db = bcr.db.BuiltinPDBDatabase()
csdb = bcr.db.BuiltinCsDatabase()
#csdb = bcr.canonical_structures.RandomForestCsDatabase()


@app.route('/')
def index():
    return render_template('index.html')


@app.route('/model', methods=['POST'])
def model():
    """Request a pdbmodel based on:

    alpha sequence
    beta sequence

    optional:
    template pdbnames
    """

    alpha = request.json.get('alphaChain')
    beta = request.json.get('betaChain')

    igc = bcr.IgComplex(alpha, beta, template_db, pdb_db)
    igc.hmmsearch(*hmms)
    igc.canonical_structures(csdb)
    igc.find_templates()#id_max=0.98)
    igc.build_structure(scwrl=False)

    pdbf = StringIO.StringIO()
    igc.save(pdbf)

    pdb = pdbf.getvalue()
    pdbf.close()

    json_data = dict()
    json_data['pdb'] = pdb

    template_sort      = sorted([chain.chain_type for chain in igc])
    json_data['templates'] = {
        'headers':   [igc[c].description for c in template_sort],
        'framework': [getattr(igc[c].templates.get('framework'), 'pdbname', None) for c in template_sort],
        'packing':   None,
        'cdrs':      [[] for c in range(3)],
    }

    for c in template_sort:
        i = 0
        for cdr, entry in sorted(igc[c].templates.items()):
            if cdr.lower() not in ('framework', 'packing'):
                json_data['templates']['cdrs'][i].append(entry.pdbname)
                i += 1

        if igc[c].templates.get('packing'):
            json_data['templates']['packing'] = igc[c].templates['packing'].pdbname

    json_data['alignment'] = {c.description: c.show_alignment(html_class='text-info') for c in igc}

    return jsonify(json_data)


@app.route('/templates')
def templates():
    """Get a list of possible templates.. Do this later."""
    return "Hello World!"


if __name__ == '__main__':
    app.run(host='0.0.0.0', debug=True)