#!/usr/bin/env python
# Copyright (C) 2015 University of Vienna
# All rights reserved.
# BSD license.
# Author: Ali Baharev <ali.baharev@gmail.com>
from __future__ import print_function
from cjson import decode as loads, encode as dumps
#from json import loads, dumps  # if you don't have cjson installed
from contextlib import contextmanager
from inspect import getmembers, getmodule, isfunction
from platform import system
import os
import sys
from sys import stdin, stdout, stderr # see also _setup_fifo_hack()
from traceback import format_exc

sys.dont_write_bytecode = True

# Whatever is printed to the sys.stdout by the called functions, it MUST NOT
# appear as it would ruin the communication; we redirect their sys.stdout to the
# log file.
LOGFILE_NAME = os.devnull

# This dictionary will be updated at the end of this file from rpc_api
registered_functions = {'shutdown': sys.exit}

#-------------------------------------------------------------------------------
# We read the input from the stdin in an infinite loop, and write our response
# to the stdout. Notifications (requests with no "id" in it) are not answered.  

def main():
    #stdin = open('data/jsontest.txt')         # <- for testing
    #tracefile = open('/tmp/json_io.log', 'w') # <- for logging
    #_setup_fifo_hack()                        # <- for interactive debugging
    #print('JSON server is up and running!')
    while True:
        request = stdin.readline().rstrip()
        #tracefile.write("Received: '{}'\n".format(request))
        with sys_stdout_redirected(LOGFILE_NAME):
            response = process(request)
        if response is not None: # None: if the request was a notification
            print(response, file=stdout)
            stdout.flush()
        #tracefile.write('{}\n'.format(response if response is not None else '(notification)'))

NOTIFICATION = 'NOTIFICATION' # Sentinel if no "id" is given in the request

def process(request_str):
    # parse request
    try:
        request = loads(request_str)
    except ValueError as valerr:
        return parse_error(str(valerr), request_str)
    except:
        return parse_error(format_exc(), request_str)
    # check request type
    if not isinstance(request, dict):
        return invalid_request('Type of request is {}'.format(type(request).__name__))
    # check protocol
    if not request.get('jsonrpc') == '2.0':
        return invalid_request('Invalid protocol')
    # get contents
    rpcid = request.get('id', NOTIFICATION)    
    method_str = request.get('method')
    method = registered_functions.get(method_str)
    if not method:
        return method_not_found('Method "{}" not found'.format(method_str), rpcid)
    params = request.get('params', [])
    # check params: only positional arguments are supported 
    if not isinstance(params, list):
        return invalid_params('If params is given, it must be a list.', rpcid)
    # call the method
    try:
        result = method(*params)
    except TypeError as ex:
        return invalid_params('Invalid parameters: {}'.format(ex), rpcid)
    except Exception:
        return internal_error(format_exc(), rpcid)
    # pack the result
    try:
        return pack(result, rpcid)
    except:
        return internal_error(format_exc(), rpcid)

#-------------------------------------------------------------------------------

def parse_error(exception_msg, request_str):
    msg = "Parse error: {}\nRequest was: '{}'".format(exception_msg, request_str)
    return _failed(msg, -32700, 'Parse error')

def invalid_request(msg):
    return _failed(msg, -32600, 'Invalid Request')

def method_not_found(msg, rpcid):
    return _failed(msg, -32601, msg, rpcid)

def invalid_params(msg, rpcid):
    return _failed(msg, -32602, msg, rpcid)

def internal_error(msg, rpcid):
    return _failed(msg, -32603, 'Internal error', rpcid)

def _failed(msg_to_log, error_code, msg, rpcid=None):
    print(msg_to_log, file=stderr)
    if rpcid is not NOTIFICATION:
        return dumps({'jsonrpc': '2.0', 
                      'error': {'code': error_code, 'message': msg}, 'id': rpcid})

def pack(result, rpcid):
    if rpcid is not NOTIFICATION:    
        return dumps({'jsonrpc': '2.0', 'result': result, 'id': rpcid})

# FIXME: What if any of the dumps calls fails? The pack is properly handled in a
# try catch, that is OK. However, if _failed raises an exception, that leaves
# the caller blocked forever... 

#-------------------------------------------------------------------------------
# We append the dir in which this script or the corresponding py2exe executable 
# lives to the sys.path. Then, we can safely import stuff from the outside 
# world such as the Python module with the RPC API. We also need it to get the 
# zip file location on Windows.
try:
    approot = os.path.dirname(os.path.abspath(__file__))
except NameError:  # We are in the main py2exe executable
    # See also http://www.py2exe.org/index.cgi/Py2exeEnvironment
    approot = os.path.dirname(os.path.abspath(sys.argv[0]))

sys.path.append(approot) # To import rpc_api on any platform, but not done yet!

if system()=='Windows':
    # Put the zip file with the Python bytecode also on the sys.path,
    # required by rpc_api. On Linux, the caller of the present script 
    # puts the developer sources on the PYTHONPATH instead.
    sys.path.append(approot + os.path.sep + 'ordering.zip')    

# TODO The ordering.zip and the hard-coded rpc_api are not generic enough

import rpc_api
MODUL = rpc_api

def to_register(o):
    return isfunction(o) and getmodule(o)==MODUL and not o.__name__.startswith('_')

registered_functions.update((name,f) for name, f in getmembers(MODUL, to_register))

#-------------------------------------------------------------------------------

@contextmanager
def sys_stdout_redirected(logfile_name):
    with open(logfile_name, 'w') as log:
        old_stdout = sys.stdout # Must keep the sys. prefix!
        sys.stdout = log
        try:
            yield
        finally:
            sys.stdout = old_stdout

#-------------------------------------------------------------------------------

def _setup_fifo_hack():
    from utils import remove_if_exists
    fifos = ('/tmp/in', '/tmp/out', '/tmp/err')
    for fifo in fifos:
        remove_if_exists(fifo)
        os.mkfifo(fifo)
    global stdin, stdout, stderr
    # Weird blocking: http://stackoverflow.com/q/5782279/341970 
    stdin  = open('/tmp/in',  'r+') # but we will only read
    stdout = open('/tmp/out', 'r+') # but we will only write
    stderr = open('/tmp/err', 'r+') # but we will only write


if __name__ == '__main__':
    main()
