#! /usr/bin/env python

import os
import uuid
import re
import subprocess
import argparse
import json
from subprocess import call
from flask import Flask, send_file, flash, send_from_directory, request, redirect, url_for, jsonify
from flask_cors import CORS
from werkzeug.utils import secure_filename

app = Flask(__name__)
CORS(app)
SAGEWS = os.path.dirname(os.path.abspath(__file__))

app.config['SAGE'] = os.path.join(SAGEWS, "..")
app.config['UPLOAD_FOLDER'] = os.path.join(app.config['SAGE'], "data")
app.config['MAX_CONTENT_LENGTH'] = 8 * 1024 * 1024   #maximum of 8MB
app.config['BASEURL'] = '/sage'
app.static_folder = app.static_url_path = os.path.join(SAGEWS, "../client/static")

def allowed_file(filename):
    return '.' in filename and filename.rsplit('.', 1)[1].lower() in set(['abi','ab1','ab!','ab', 'json', 'fa'])

@app.route('/upload', methods=['POST'])
def generate():
    uuidstr = str(uuid.uuid4())

    # Get subfolder
    sf = os.path.join(app.config['UPLOAD_FOLDER'], uuidstr[0:2])
    if not os.path.exists(sf):
        os.makedirs(sf)

    if request.form['sample'] == 'sample':
        fexpname = os.path.join(SAGEWS, "sample.abi")
        genome = os.path.join(SAGEWS, "sample.fa")
    else:
        # Experiment
        if 'experiment' not in request.files:
            return jsonify(errors = [{"title": "Experiment file is missing!"}]), 400
        fexp = request.files['experiment']
        if fexp.filename == '':
            return jsonify(errors = [{"title": "Experiment file is missing!"}]), 400
        if not allowed_file(fexp.filename):
            return jsonify(errors = [{"title": "Experiment file has incorrect file type!"}]), 400
        fexpname = os.path.join(sf, "sage_" + uuidstr + "_" + secure_filename(fexp.filename))
        fexp.save(fexpname)

        # Genome
        if request.form['refType'] == 'genome':
            genome = request.form['genome']
            if genome == '':
                return jsonify(errors = [{"title": "Genome index is missing!"}]), 400
            genome = os.path.join(app.config['SAGE'], "fm", genome)
        elif request.form['refType'] == 'fasta':
            if 'fasta' not in request.files:
                return jsonify(errors = [{"title": "Fasta file is missing!"}]), 400
            fafile = request.files['fasta']
            if fafile.filename == '':
                return jsonify(errors = [{"title": "Fasta file is missing!"}]), 400
            if not allowed_file(fafile.filename):
                return jsonify(errors = [{"title": "Fasta file has incorrect file type!"}]), 400
            genome = os.path.join(sf, "sage_" + uuidstr + "_" + secure_filename(fafile.filename))
            fafile.save(genome)
        elif request.form['refType'] == 'trace':
            if 'trace' not in request.files:
                return jsonify(errors = [{"title": "Wildtype Chromatogram file is missing!"}]), 400
            tracefile = request.files['trace']
            if tracefile.filename == '':
                return jsonify(errors = [{"title": "Wildtype Chromatogram file is missing!"}]), 400
            if not allowed_file(tracefile.filename):
                return jsonify(errors = [{"title": "Wildtype Chromatogram file has incorrect file type!"}]), 400
            genome = os.path.join(sf, "sage_" + uuidstr + "_" + secure_filename(tracefile.filename))
            tracefile.save(genome)
        else:
            return jsonify(errors = [{"title": "No input reference provided!"}]), 400

    # Run sage
    outfile = os.path.join(sf, "sage_" + uuidstr + ".json")
    alnfile = os.path.join(sf, "sage_" + uuidstr + ".align")
    logfile = os.path.join(sf, "sage_" + uuidstr + ".log")
    errfile = os.path.join(sf, "sage_" + uuidstr + ".err")
    with open(logfile, "w") as log:
        with open(errfile, "w") as err:
            blexe = os.path.join(app.config['SAGE'], "./src/sage")
            return_code = call([blexe,'-g', genome,'-o',outfile,'-a',alnfile,fexpname], stdout=log, stderr=err)
    if return_code != 0:
        errInfo = "!"
        with open(errfile, "r") as err:
            errInfo = ": " + err.read()
        return jsonify(errors = [{"title": "Error in running sage" + errInfo}]), 400
    return jsonify(data = json.loads(open(outfile).read()))


@app.route('/')
def root():
    return send_from_directory(os.path.join(SAGEWS, "../client"),"index.html"), 200


if __name__ == '__main__':
    app.run(host = '0.0.0.0', port=3300, debug = True, threaded=True)

                                                                        
