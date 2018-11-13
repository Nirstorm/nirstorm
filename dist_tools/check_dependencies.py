"""
Check matlab backward compatibility for all functions used in the package.

Parse all matlab functions used in all .m files in the nirstorm directory 
(recursive) and query matlab reference document to know when the function 
was released. 
"""
import os
import os.path as op
import re
import requests
import subprocess
from collections import namedtuple, defaultdict

DEBUG = True

func_re = re.compile('.*?([a-zA-Z]\w+)\([^)]*\)[^.,].*?')
release_re = re.compile('Introduced (in|before) (\w+) *<')

print 'Retrieving matlab reference documentation path...'
MATLAB_DOC_PATH_CMD = ['matlab', '-nosplash', '-nodesktop', '-noFigureWindows',
                       '-r' , 'try; fprintf("matlabroot:%s@", matlabroot); '\
                              'catch; display(\'Error\'); exit(1); end; exit(0);']
print ' '.join(MATLAB_DOC_PATH_CMD)
MATLAB_DOC_PATH_RE = re.compile('matlabroot:(.*)@.*$')
if 0:
    mat_output = subprocess.check_output(MATLAB_DOC_PATH_CMD)
    # print 'matlab output:'
    # print mat_output
    doc_re_result = MATLAB_DOC_PATH_RE.findall(mat_output)
    if len(doc_re_result) == 0:
        raise Exception('Could not retrieve matlab root path')
    if len(doc_re_result) > 1:
        raise Exception('Error, got multiple matlab root paths:\n' + \
                        '\n'.join(doc_re_result))
    
    MATLAB_ROOT = doc_re_result[0]
    if not op.exists(MATLAB_ROOT):
        raise Exception('Retrieved matlab root path does not exist: ' + MATLAB_ROOT)
    MATLAB_HELP_PATH = op.join(MATLAB_ROOT, 'help')
else:
    MATLAB_HELP_PATH = '/home/tom/Projects/Software/matlab_official_2017a/help'
    
print 'Parsing matlab doc path ',  MATLAB_HELP_PATH, '...'

matlab_help_fns = defaultdict(list)
for root, dirs, fns in os.walk(MATLAB_HELP_PATH):
    for fn in fns:
        if fn.endswith('.html'):
            matlab_help_fns[op.splitext(op.basename(fn))[0]].append(op.join(root, fn))


LocalFileResponse = namedtuple('LocalFileResponse', 'ok text toolbox')

def get_matlab_help_doc(func_name):

    if matlab_help_fns.has_key(func_name):
        doc_fns = matlab_help_fns[func_name]
        if len(doc_fns) > 1:
            print 'Warning, more than one version of ', func_name, 'found:'
            print '\n'.join(['  - %s' %d for d in doc_fns])

        # Pick core version if exists: 
        doc_fn = None
        for dfn in doc_fns:
            if op.join('matlab', 'ref') in dfn:
                doc_fn = dfn
                toolbox = 'core'
                break
        if doc_fn is None:
            doc_fn = doc_fns[0]
            toolbox = op.relpath(doc_fn, MATLAB_HELP_PATH).split(op.sep)[0]
                
        print 'Using local doc file:', doc_fn, '(tbx:', toolbox ,')'
        with open(doc_fn) as fdoc:
            return LocalFileResponse(True, fdoc.read(), toolbox)
    else:
        return LocalFileResponse(False, '', None)
        # To get doc from mathworks website:
        # url = matlab_doc_url_pat % func_name
        # print 'Trying URL:', url
        # return requests.get(url, timeout=1)

def describe_matlab_func(func_name, searched):
    searched.add(func_name)
    date_tag = None
    toolbox = None
    if DEBUG:
        print 'Checking doc of %s ...' % func_name
    try:
        resp = get_matlab_help_doc(func_name)
        toolbox = resp.toolbox
        if resp.ok:
            # Extract release date
            rr = release_re.search(resp.text)
            if rr is not None:
                if rr.group(1) == 'in':
                    date_tag = rr.group(2)
                else: #before
                    date_tag = 'b4_' + rr.group(2)
            else:
                if DEBUG:
                    print func_name, ': rdate not found in doc'
        else:
            if DEBUG:
                print func_name, ': func doc not found'
    except requests.exceptions.ConnectionError:
        if DEBUG:
            print func_name, ': func doc not found online (connection error)'
    except requests.exceptions.ReadTimeout:
        if DEBUG:
            print func_name, ': func doc not found online (time out)'
    
    return date_tag, toolbox

def can_be_mat_func(func_name):
    """
    Return False if given function name cannot be a core matlab function.
    Eg:
      - starts with process_ or bst_ (brainstorm-specific)
      - starts with nst_ (nirstorm-specific)
    """
    return not (func_name.startswith('process_') or \
                func_name.startswith('bst_') or \
                func_name.startswith('nst_') or \
                func_name.startswith('mfip_'))
    
from glob import glob
rdates = defaultdict(set)
toolboxes = defaultdict(set)
ignore = set()

for mat_fn in glob(op.join(op.dirname(op.realpath(__file__)), '../**/*.m')):
    print 'Parsing %s...' % mat_fn
    with open(mat_fn) as fmat:
        func_names = func_re.findall(fmat.read())
        func_info = [(f,describe_matlab_func(f,ignore)) for f in func_names
                   if f not in ignore and can_be_mat_func(f)]
        for f,(rdate, tbx) in func_info:
            if rdate is not None:
                rdates[rdate].add(f)
            if tbx is not None:
                toolboxes[tbx].add(f)

print '\n-- Functions sorted by release date --\n'
for rdate, funcs in rdates.items():
    print rdate
    print ', '.join(sorted(funcs))
    print ''

print '\n-- Functions sorted by toolbox --\n'
for tbx, funcs in toolboxes.items():
    print tbx
    print ', '.join(sorted(funcs))
    print ''
    
