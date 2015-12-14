#!/usr/bin/env python
#
# Copyright (C) 2014, 2015 University of Vienna
# All rights reserved.
# BSD license.
# Author: Ali Baharev <ali.baharev@gmail.com>
from __future__ import print_function
from gzip import GzipFile
from xml.etree import ElementTree as ET
import os
from os.path import isdir
import shutil
import subprocess
import sys
from tempfile import gettempdir
from zipfile import ZipFile
from utils import DATADIR

JM_PYTHON  = '/home/ali/jmodelica_install/bin/jm_python.sh'
FLATTEN_PY = '/home/ali/ws-pydev/structure-reconstruction/fmux_creator.py'
# Directories must have the trailing os.sep ('/' or '\')
TMPDIR     = gettempdir()+os.sep+'FMUX'+os.sep

def flatten(model_file, class_name):
    create_fmux(model_file, class_name)
    model_name = model_file.rpartition('.')[0] # works for both .mo and .mop
    unpack_fmux(model_name, class_name)
    clean_tmpdir(model_name)

def create_fmux(model_file, class_name):
    shutil.rmtree(TMPDIR, ignore_errors=True)
    os.mkdir(TMPDIR)
    shutil.copy(DATADIR+model_file, TMPDIR)
    cmd = [JM_PYTHON, FLATTEN_PY, model_file, class_name]
    if subprocess.call(cmd, cwd=TMPDIR):
        sys.exit(1)  

def unpack_fmux(model_name, class_name):
    mangled_name = class_name.replace('.', '_')
    with ZipFile(TMPDIR+mangled_name+'.fmux') as f:
        f.extractall(path=TMPDIR)
    # The get_ampl_name relies on these extensions (.xml and .xml.gz)
    flatmodel_xml = TMPDIR+model_name+'.xml'
    os.rename(TMPDIR+'modelDescription.xml', flatmodel_xml)
    compress(flatmodel_xml)
    # Move the flattened .mof files from the ./sources/ to the temp dir
    SRCDIR = TMPDIR+'sources'+os.sep
    shutil.move(SRCDIR+class_name+'.mof', TMPDIR+model_name+'.mof' )
    suffix = '_transformed.mof'
    shutil.move(SRCDIR+class_name+suffix, TMPDIR+model_name+suffix )

def compress(filename):
    with open(filename, 'rb') as orig, GzipFile(filename+'.gz', 'wb') as out:
        out.writelines(orig)

def clean_tmpdir(model_name):
    to_rm=[TMPDIR+e for e in os.listdir(TMPDIR) if not e.startswith(model_name)]
    for entry in to_rm:
        os.remove(entry) if os.path.isfile(entry) else shutil.rmtree(entry)

def get_ampl_name(xml_filename):
    ampl_name = TMPDIR+get_model_name(xml_filename)+'.ampl'
    if not isdir(TMPDIR): 
        os.mkdir(TMPDIR) # or open(ampl_name) throws otherwise 
    return ampl_name 

def get_model_name(xml_filename):
    # A hideous way to get the model name
    xml = os.path.basename(xml_filename)
    assert xml.endswith('.xml') or xml.endswith('.xml.gz'), xml
    crop = len('.xml') if xml.endswith('.xml') else len('.xml.gz')
    return xml[:-crop]

def get_etree_without_namespace(xml_file):
    xml_file = DATADIR + xml_file
    if xml_file.endswith('.gz'):
        xml_file = GzipFile(xml_file)
    tree = ET.parse(xml_file)
    # Hackish removal of the namespace
    for elem in tree.iter():
        i = elem.tag.find('}')
        elem.tag = elem.tag[i+1:]
        #print(elem.tag, '  ', elem.attrib)
    return tree

if __name__=='__main__':
    if len(sys.argv)!=3:
        print('Usage: %s model_file class_name' % sys.argv[0], file=sys.stderr)
        sys.exit(0)
    flatten(sys.argv[1], sys.argv[2])
