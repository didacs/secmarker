#!/usr/bin/env python
import sys,os,Secmarker
test_folder = os.path.dirname(os.path.realpath(__file__))+'/test'

print
print '** SECMARKER test **'
print '********************'
print
print 'run: ./Secmarker.py -t '+ test_folder +'/test.fa -o '+ test_folder +'/output'

args = Secmarker.parse_opts( '-t '+ test_folder +'/test.fa -o '+ test_folder +'/output')

if Secmarker.main(args):
    print 'cmsearch test: OK'
else:
    print 'cmsearch test: ERROR!'
    print 'Something failed in the cmsearch phase test...'
    print 'exiting'
    sys.exit()

print 'run: ./Secmarker.py -t '+ test_folder +'/test.fa -o '+ test_folder +'/output -plot'
args = Secmarker.parse_opts( '-t '+ test_folder +'/test.fa -o '+ test_folder +'/output -plot')
if Secmarker.main(args):
    print 'RNAplot test: OK'
else:
    print 'WARNING: Secmarker did not work with -plot option. Install ViennaRNA package to use it...'
