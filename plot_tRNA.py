import sys
from MMlib import bash
from Secmarker import rnaplot_arm_colors,cca_coords,expand_boundingbox

"""usage: python plot_tRNA.py input
input is the trnasec.ss file
the idea of this script is to be able to modify the .ss file in order to generate a plot from it

the tRNA sequence and structure is obtained from the ss file

rf  AAAAAAAA,,DDDD.DD====DDDDDDC.CCCC==ACD==CCCC...CVVVVV.VVV====VVV.VVVVV,TTTTT=======TTTTTAAAAAAAAdcca
seq UGGGAGGUUAGUGUACCUGGUGGGCACCACAGGCUUCAACCCUGAUUGACGGUUUUGAUGACAACGCCGUGGGAGGUUCGAUUCCUUCACCUUCUCGCCA
ss  .(((((((..((((.((....))))))(.((((.......))))...)(((((.(((....))).))))).(((((.......)))))))))))).....

all pairs in the ss must have any of the ADCVT letter in the rf line, each pair will be colored according to the arm

options in the rf line:
the anticodon position is marked by ACD (uppercase) in the rf line (note that Secmarker output is in lowercase)
the discr base position is marked by d
if the rf line contains cca, the plot will contain the CCA

"""

trna_ss = sys.argv[1]

# functions, classes and variables
class symmetrical_hash(dict):
    """ ........... """
    def get(self, k1, k2):
        try:
            return self[k1][k2]
        except:
            return self[k2][k1]

category_per_pair=symmetrical_hash()
category_per_pair['U']={'A':'a'}
category_per_pair['T']={'A':'a'}
category_per_pair['G']={'A':'k', 'C':'c', 'U':'w', 'T':'w'}
category_per_pair_hierarchy='cawko'

def pairs(ss_string):
    """ Given a string with secondary structure  -- in which pairs are represented by any parenthesis ([{<  -- returns a list of tuples with 0 based positions of each pair, from 5' to 3'"""
    height = 0;     stem = {};     pairs = []; #    pos_rf=0;     pos_t=0                                                                                                                                                                    
    for pos, ss_char in enumerate(ss_string):
        if ss_char in '([{<':
            stem[ height ] = pos
            height += 1
        elif ss_char in ')]}>' and height > 0:
            height -= 1
            if height < 0: break
            paired_pos = stem[height]
            pairs.append( (paired_pos, pos) )
    return sorted(pairs, key = lambda x: x[0])

def target_ss_pairs(rf,ss,seq):
    target_ss_pairs = []
    for i,j in pairs(ss):
        if not rf[i] == rf[j]:
            raise Exception, 'ERROR: wrong hand_rf annotation: %s-%s for this pair: %s-%s' % (rf[i], rf[j], i,j)
        nt_i = seq[i]
        nt_j = seq[j]
        try:
            category=category_per_pair.get(nt_i,nt_j)
        except KeyError:
            category='o'
        target_ss_pairs.append( (i+1,j+1,category,rf[i]) )
    return target_ss_pairs


d = {}
for i in open(trna_ss):
    if not i.rstrip().split(): continue
    if i.startswith('>>'):
        trna_id = int(i.rstrip().split(':')[-1])
        if i in d:
            raise Exception, 'ERROR: tRNA identifiers not unique'
    d.setdefault(trna_id, []).append(i.strip())

for trna_id,v in d.items():

    rf  = v[-3]
    seq = v[-2].replace('T','U')
    ss  = v[-1]

    # control paired parentheses                                                                                                                                                                                                             
    if ss.count('(') != ss.count(')'):
        raise Exception, 'ERROR: malformated ss, unbalanced parentheses in '+ss

    pre = '--pre "/dbasemark { newpath 1 sub coor exch get aload pop fsize 1.75 div 0 360 arc stroke} bind def /outlinecolor {0.6 setgray} bind def /paircolor {0.6 setgray} bind def '
    for i,j,category,arm in target_ss_pairs(rf,ss,seq):
        sat=.5
        if not category in 'ko':
            if category ==  'w': sat=.25
            pre +=  ' %s %s %s %s colorpair' % (i,j,rnaplot_arm_colors[arm],sat)
    # anticodon position highlight in pre string
    anticodon_pos = rf.find('ACD') + 1
    if anticodon_pos > -1:
        for i in range(3):
            pre += ' %s cmark' % (anticodon_pos + i)
    # cmark in discriminator base
    dbase = rf.find('d') + 1
    if dbase > -1:
        pre += ' %s dbasemark' % (dbase)
    pre += '"'


    rnaplot_input =  '>tRNA_%s\n%s\n%s' % (trna_id,seq,ss)
    b = bash('echo "'+ rnaplot_input +'" | RNAplot -t 1 '+ pre )
    if 'ERROR' in b[1]: raise Exception, b[1]
    ps_file = 'tRNA_%s_ss.ps' % trna_id
    # modify coords of CCA 3' end, if 'cca' is defined in rf string
    if rf.find('cca') > -1:
        cca_coords(ps_file)
    # expand boundingbox. This solves the problem of having the image cut at the cmarks in the anticodon                                                                                                                                     
    expand_boundingbox(ps_file)
    png_file = 'tRNA_%s_ss.png' % trna_id
    b = bash('convert -density 150 '+ ps_file +' '+ png_file)
    if b[0]: raise Exception, b[1]
