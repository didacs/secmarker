#!/usr/bin/env python
import os,sys
import argparse

setup_dir = os.path.dirname(os.path.realpath(__file__))
library = setup_dir +'/library'
bin_folder = setup_dir +'/bin'

sys.path.insert(0, library)
from MMlib import bash,bbash

parser = argparse.ArgumentParser(description="setup.py checks if cmsearch, esl-sfetch and RNAplot programs are available. It also checks the version of Infernal (v1.1.1 required). It then creates soft links inside the ./bin folder. Without arguments cmsearch, esl-sfetch and RNAplot should be in your $PATH, otherwise the path to the executables can be provided with the corresponding arguments.")
parser.add_argument('--cmsearch',        help='path to cmsearch executable: required INFERNAL 1.1.1')
parser.add_argument('--sfetch',        help='path to esl-sfetch executable (from Infernal package)')
parser.add_argument('--RNAplot',        help='path to RNAplot executable: recomended RNAplot >= 2.1.2')
args = parser.parse_args()

def check_infernal_version(b):
    if b[0]:
        if 'command not found' in b[1]:
            raise Exception, 'ERROR: cmsearch program not available'
        else:
            raise Exception, b[1]
    for i in b[1].split('\n'):
        if i.startswith('# INFERNAL'):
            if not '# INFERNAL 1.1.1' in i:
                raise Exception, 'ERROR: Infernal version must be 1.1.1, current installed: '+ i
            return True

print 'CHECK cmsearch...'
# create cmsearch link into bin folder
if args.cmsearch:
    cmsearch = args.cmsearch
else:
    b = bash('which cmsearch')
    if b[0] or not b[1]: raise Exception, 'ERROR: cmsearch is not in your $PATH, please use "--cmsearch /path/to/cmsearch"'
    cmsearch = b[1]

b = bash(cmsearch +' -h')
if check_infernal_version(b):
    bbash('ln -sf '+ cmsearch +' '+bin_folder)
    print 'link created: ./bin/cmsearch ->', cmsearch


print 'CHECK esl-sfetch...'
# create --esl-sfetch link into bin folder
if args.sfetch:
    sfetch = args.sfetch
else:
    b = bash('which esl-sfetch')
    if b[0] or not b[1]: raise Exception, 'ERROR: esl-sfetch is not in your $PATH, please use "--sfetch /path/to/esl-sfetch". \
esl-sfetch is part of the infernal package, it should be inside the binaries/ folder, or the src/ folder."'
    sfetch = b[1]

bbash('ln -sf '+ sfetch +' '+bin_folder)
print 'link created: ./bin/esl-sfetch ->', sfetch


print 'CHECK RNAplot...'    
# create cmsearch link into bin folder                                
if args.RNAplot:
    rnaplot = args.RNAplot
else:
    b = bash('which RNAplot')
    if b[0] or not b[1]: raise Exception, 'ERROR: RNAplot is not in your $PATH, please provide the /full/path/to/RNAplot as --RNAplot option'
    rnaplot = b[1]
bbash(rnaplot +' -h')
bbash('ln -sf '+ rnaplot +' '+bin_folder)
print 'link created: ./bin/RNAplot ->', rnaplot

print 'DONE!'
