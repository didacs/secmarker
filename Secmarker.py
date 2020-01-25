#!/usr/bin/env python
# AUTHOR: Didac Santesmasses <didac.santesmasses@crg.eu>

##########################################################################
#
#    This file is part of Secmarker. http://secmarker.crg.eu/
#
#    Secmarker is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#
#    Secmarker is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#
#    You should have received a copy of the GNU General Public License
#    along with Secmarker.  If not, see <http://www.gnu.org/licenses/>.
#
##########################################################################

__VERSION__="secmarker-0.4a" 

import sys,os
import argparse
sys.path.insert(0, os.path.dirname(os.path.realpath(__file__)) +'/library')
from MMlib import *
os.environ['PATH'] = os.path.dirname(os.path.realpath(__file__)) +'/bin:'+os.environ['PATH']

def get_parser():
    """Read the arguments for the program"""
    parser = argparse.ArgumentParser(description='search for selenocysteine tRNA (tRNA-Sec). The presence of tRNA-Sec in a genome indicates the use of selenocysteine by the organism')
    parser.add_argument('-t','--target', dest='target', help='target nucleotide sequence', required = True)
    parser.add_argument('-o', help='output folder', default='.')
    parser.add_argument('-plot', help='generate the tRNA plot(s). For each tRNA prediction, two files (.png and .ps) will be generated in the folder images/', action="store_true")
    parser.add_argument('-T','--infernal_score', dest='T', help='score threshold used by cmsearch (Infernal score)', default=40, type=float)
    parser.add_argument('-F','--filtering_score', dest='F', help='score threshold used to filter tRNAs. If a tRNA hit has the UCA anticodon or a G73, -T is used for filtering', default=55, type=int)

    parser.add_argument('-cpu', help='number of available CPUs [default 1]', default='1')
    parser.add_argument('-cca', help='include the 3 residues downstream the discriminator base if available (useful to see whether the CCA sequence is genome encoded)', action="store_true")
    parser.add_argument('-AT12', help='force a 12 bp AT-stem, that is a 7/5 fold tRNA', action="store_true")
    parser.add_argument('-O','--options',dest='O', help='additional options for cmsearch', default='', type=str)
    return parser


def parse_opts(x):
    """To be used when the program is imported as a module:
import Secmarker
args = Secmarker.parse_opts('-y your -o options -f for -t the -p program')
Secmarker.main(args)
"""
    parser = get_parser()
    args = parser.parse_args(x.split())
    return args


def main(args):
########## PROGRAM START
    # output folder 
    args.o = Folder( args.o.rstrip('/')+'/' )
    # library folder
    global secmarker_library; secmarker_library = os.path.dirname(os.path.realpath(__file__))+'/library/'
    if not os.path.exists(secmarker_library):
        raise Exception, "ERROR folder %s does not exist." % secmarker_library

    # cm file
    cm = secmarker_library+'tRNAsec_db.cm'
    if not os.path.isfile(cm):
        raise Exception, "ERROR cm:%s file not found." % cm
    # check if the required programs are available: cmsearch, esl-sfetch and RNAplot
    check_depencies(args)

    # check target file presence
    if not os.path.isfile(args.target):
        raise Exception, "ERROR target:%s file not found." % args.target

    # generate a soft link of target file
    # this will solve the problem of not having permission in the target folder when writting the esl-sfetch index
    target_link = args.o +'genome_link.fa'
    if not os.path.isfile(target_link):
        b = bash('ln -s '+ os.path.abspath(args.target) +' '+ target_link)
        if b[0]: raise Exception, b[1]
    if not os.path.exists(os.readlink(target_link)):
        raise Exception, 'ERROR the soft link for the target sequence in the output folder:%s was not correctly created...' % (args.o)
    # keep original target path and use the soft link as the main target file
    user_input_target = args.target
    args.target = target_link

    # generate esl-sfetch .ssi index file
    if not os.path.isfile(args.target +'.ssi'):
        b=bash('esl-sfetch --index '+ args.target)
        if b[0]: raise Exception, b[1]
        if not os.path.isfile(args.target +'.ssi'):
            raise Exception, "ERROR unable to generate target index .ssi file:%s ." % args.target+'.ssi'


    ########## computing chromosome lengths if necessary, loading them
#    chromosome_length_file = args.o+'target.length'
#    if not is_file(chromosome_length_file):
#        bbash('fastalength '+ args.target + ' > '+ args.o +'chrom_lengths && mv '+ args.o +'chrom_lengths '+ chromosome_length_file)

#    max_length=63
#    load_chromosome_lengths(chromosome_length_file)#,  max_length) # this will make the chromosome lengths available to selenoprofiles module and to MMlib; it's also checking the maximum number of characters for the fasta headers
#    global chromosome_lengths; chromosome_lengths = get_MMlib_var('chromosome_lengths')
    global chromosome_lengths; chromosome_lengths = {}

    ########## run infernal
    cmsearch_tmp = args.o + 'infernal.out.tmp'
    cmsearch_outfile = args.o + 'infernal.out'
    if os.path.isfile(cmsearch_outfile):
        # if infernal file exists, check if the target is the same
        target_is_the_same, infernal_target = check_target(cmsearch_outfile, args.target)
        if not target_is_the_same:
            raise Exception, 'ERROR: output folder contains infernal outputs from another target: %s' % infernal_target
            
    else:
        b = bash('cmsearch --cpu '+ args.cpu +' --notextw -T '+ str(args.T) +' -o '+ cmsearch_tmp +' '+ args.O +' '+ cm +' '+ args.target)
        if b[0]: raise Exception, b[1]
        bbash('mv '+ cmsearch_tmp +' '+cmsearch_outfile)
        
    
    # parse infernal output
    hits_list = []
    infernal = cmsearch_outfile


    # remove overlapping hits from the different models to produce a non-redundant list of hits
    # originally this was called after parsing and filtering the hits, but it produced a bug when two adjacent hits were found. 
    # parsed hits are extended by one base to get the discriminator base, so after extending the first of two adjacent hits, it overlapped the 2nd one
    infernal_hits_list = parse_infernal(infernal).all()

    # for each hit generate a gene object corresponding to the anticodon,
    # the hit.id is a stable unique identifier shared by the hit and the corresponding acd gene, here it is stored in a dict to be used after merging acd genes and retrieve back the corresponding hits
    id_dict = {}
    acd_list = []
    for hit in infernal_hits_list:
        acd = hit.hand_rf.find('acd')
        acd_in_seq = hit.alignment.position_in_seq('t', acd+1) # 1 based

        if hit.strand == '+':
            left = acd_in_seq - 1
            acd_gene = hit.extend(-left) # copy of hit, keeping its attributes but modifying the coords
            acd_gene.exons[0][1] = acd_gene.exons[0][0] + 2
        else:
            right = acd_in_seq - 1
            acd_gene = hit.extend(0,-right) # copy of hit, keeping its attributes but modifying the coords
            acd_gene.exons[0][0] = acd_gene.exons[0][1] - 2

        if hit.id in id_dict: raise Exception, 'ERROR this hit.id:%s is not unique!' % (hit.id)
        id_dict[hit.id] = hit
        acd_list.append(acd_gene)

    # remove redundancy taking into account only the anticodon
    def mode_function(a,b):
        return sorted([a, b], key=lambda x: x.score, reverse=True)[0]
    nr_acd_list = merge_genes(acd_list, strand=True, phase=False, mode = mode_function)
    nr_hits_list = [id_dict[x.id] for x in nr_acd_list]

#    nr_hits_list = merge_genes(infernal_hits_list, strand=True, phase=False, mode = mode_function)

    # parse non-redundant hits
    hit_index = 0
    hits_list = [] # final list of accepted hits
    for hit in nr_hits_list:
        hit_index += 1

        # define some hit attributes
        hit.indexx = hit_index # merge_genes function called below deletes the .index attribute
        hit.target = args.target
        hit.cmali = alignment()
        infernal_stk = secmarker_library + hit.query.chromosome +'.stk'
        hit.cmali.load_stockholm( infernal_stk )
        hit.cmali.add('rf', hit.cmali.rf.replace('.','-'))
        assert hit.cmali.check_length()
        hit.target_ss_pairs = target_ss_pairs(hit) # pairs target sequence based, not hit based

        # refine structure (adding additional good pairs in loops). This function requires hit.target_ss_pairs, and updates it.
        refine_structure(hit)

        # get the discriminator base (position 73)
        # this step also includes the extension to get the cca sequence, if -cca option is active
        hit = hit.discriminator_base(args, with_cca=args.cca)

        # remove "X" residues from the infernal output
#        if "X" in hit.alignment.seq_of('t').upper():
#            gene_seq = hit.sfetch().upper().replace('T','U')
#            try:
#                hit.remove_Xs(gene_seq)
#            except:
#                print hit
#                raise

        # get anticodon position
        try:
            acd = hit.hand_rf.find('acd')
        except AttributeError:
            raise Exception, 'ERROR %s generated with an obsolete cm, remove that file or select a new output folder with option "-o new_folder"' % (infernal)

        if acd < 0: continue         # skip hits that do not span the anticodon
        hit.acd = hit.alignment.seq_of('t')[ acd:acd+3 ].upper()


        # two different score cut-offs are used here, depending on the anticodon and the discr base
        # default cut-off -> -F (--filtering_score)
        # genes with the UCA anticodon or G73 -> -T (--infernal_score)
        score_cutoff = args.F
        if hit.acd == 'UCA' or hit.discriminator_base == 'G':
            score_cutoff = args.T

        # skip hits that do not pass the score cut-off
        if hit.score < score_cutoff: continue

        # hits that passed the first filtering
        hits_list.append(hit)


    # sort hits by score
    hits_list.sort(key = lambda x: x.score, reverse = True)

    #### filter tRNASec
    trnasec_list = []
    for hit in hits_list:
        if filter_tRNASec(hit):

            hit_id(hit) # add (e)ukaryotic, (b)acterial or (a)rchaeal to hit id label
            hit.target_ss, hit.target_hand_rf = target_ss(hit) # get target based structure, not hit based
            # extend target_ss and target_hand_rf attributes, both are used in draw_tRNA, and are included in the .ss file in output
            extension = len(hit.extended_seq)-len(hit.target_ss)
            hit.target_ss += '.'*extension
            hit.target_hand_rf += '.'*extension

            # produce graphical output
            if args.plot:
                draw_tRNA(hit, args.o)

            # final list of candidates
            trnasec_list.append(hit)


            # get seq A-T stem plus extension
#            if 'T' in hit.target_ss_pairs[-1][-1]:
#                last_j_T_stem = hit.target_ss_pairs[-1][1]
#                right = 50
#                left = 0
#                if hit.strand == '-':
#                    left,right = right,left
#                extended_gene = hit.extend(left,right)
##                extended_gene = hit.extend(5,5)
#
#                ext_seq = extended_gene.sfetch()
#                trna_type = hit.query.chromosome[0]+'tRNAsec'
#                sp = hit.target.split('/')[-2]
#                header = '.'.join([sp,trna_type,str(hit.indexx)])+' '+hit.header()+' infernal_score:'+ str(hit.score) +' discr_base:'+hit.discriminator_base
#                comment = 'tRNAsec:%s model:%s infernal_score:%s target:%s discr_base:%s' % (hit.indexx,hit.query.chromosome,hit.score,hit.target,hit.discriminator_base)
#                print '>%s from:%s\n%s' % (header, last_j_T_stem, ext_seq.replace('T','U')[last_j_T_stem:])

#    sys.exit() # skip output phase
    #### output
    if __name__=='__main__' and not trnasec_list:
        sys.stderr.write('Computation finished. No tRNA-Sec detected.\n')
        sys.exit()

    output_gff = args.o + 'trnasec.gff'
    output_fa = args.o + 'trnasec.fa'
    output_ss = args.o + 'trnasec.ss'
    output_gff_fh = open(output_gff,'w')
    output_fa_fh = open(output_fa,'w')
    output_ss_fh = open(output_ss,'w')

    for hit in trnasec_list:
        hit.target = user_input_target # original target, not the soft link
        if hit.target.endswith('/genome.fa'):
            sp = hit.target.split('/')[-2]
        else:
            sp = hit.target.split()[0]

        trna_type = hit.query.chromosome[0]+'tRNAsec'
        truncated = hit.ss.startswith('~') or hit.ss.endswith('~')
        header = '.'.join([sp,trna_type,str(hit.indexx)])+' '+hit.header()+' infernal_score:'+ str(hit.score) +' truncated:'+str(truncated) +' discr_base:'+hit.dbase +' anticodon:'+ hit.acd
        comment = 'tRNAsec:%s model:%s infernal_score:%s truncated:%s target:%s discr_base:%s anticodon:%s' % (hit.indexx,hit.query.chromosome,hit.score,truncated,hit.target,hit.dbase,hit.acd)
        
        print >> output_fa_fh, '>%s\n%s' % (header, hit.extended_seq)
        print >> output_gff_fh, hit.gff( tag = 'tRNA', comment = comment , program = 'secmarker')
        print >> output_ss_fh, hit_nice_output(hit)

        if __name__=='__main__':
            write( hit.gff( tag = 'tRNA', comment = comment, program = 'secmarker' ) + '\n' )

    output_gff_fh.close()
    output_fa_fh.close()

    if __name__=='__main__':
        sys.stderr.write('Computation finished. %s tRNA-Sec detected.\n' % len(trnasec_list))

    else:
        return trnasec_list


#### classes
##################################################


#### functions and variables
##################################################

# add sfetch method to infernalhit class
def sfetch(self):
    target = self.target
    if not os.path.isfile(target +'.ssi'):
        try:
            bbash('esl-sfetch --index '+ target)
        except:
            bbash('ln -sf '+ target +' '+ args.o +'genome_link.fa')
            target = args.o +'genome_link.fa'
            if not os.path.isfile(target +'.ssi'):
                bbash('esl-sfetch --index '+ target)
    sequence = ''
    ifrev = int(self.strand=='-')
    for exon_index in range(len(self.exons)):
        start, stop = self.exons[exon_index]
        b = bash('esl-sfetch '+ '-r '*ifrev +'-c %s-%s "%s" "%s"' % (start, stop, target, self.chromosome))
        if b[0]: raise Exception, b[1]
        sequence += ''.join( [x for x in b[1].split('\n') if not x.startswith('>')] )
    return sequence
setattr(infernalhit, 'sfetch', sfetch)

# add discriminator_base method to infernalhit class
def discriminator_base(self, args, with_cca=False):
    # extend gene
    left,right = define_extension(self, args, with_cca=with_cca)
    extended_gene = self.extend(left,right, inplace=False)

    # compute chromosome lengths
    if not self.chromosome in chromosome_lengths:
        b = bash('esl-sfetch %s "%s"' %( args.target, self.chromosome))
        seq = ''.join( b[1].split('\n')[1:] )
        chromosome_lengths[self.chromosome] = len(seq)

    # check if boundaries were modified, meaning the extension is out of boundaries
    if extended_gene.check_boundaries(chromosome_length=chromosome_lengths[self.chromosome]):
        has_cca = False
        has_dbase = False 
        if with_cca: # the cca was asked and it is out of boundaries, lets try just the dbase
            left,right = define_extension(self, args, with_cca=False)
            extended_gene = self.extend(left,right, inplace=False)
            if not extended_gene.check_boundaries(chromosome_length=chromosome_lengths[self.chromosome]):
                has_dbase = True # the discr base is now contained in the gene
    else:
        has_cca = with_cca
        has_dbase = True

    # define dbase, dbase_position, cca, has_dbase and has_cca attrs
    extended_seq = extended_gene.sfetch().upper().replace('T','U')
    extended_gene.extended_seq = extended_seq
    extended_gene.has_cca = has_cca
    extended_gene.has_dbase = has_dbase
    if extended_gene.has_cca:
        extended_gene.dbase = extended_seq[-4]
        extended_gene.cca =  extended_seq[-3:]
        extended_gene.dbase_position = len(self.sequence()) + max(left,right) - 3
    elif extended_gene.has_dbase:
        extended_gene.dbase =extended_seq[-1]
        extended_gene.cca = '?'
        extended_gene.dbase_position = len(self.sequence()) + max(left,right)
    else:
        extended_gene.dbase = '?'
        extended_gene.cca = '?'
    return extended_gene
setattr(infernalhit, 'discriminator_base', discriminator_base)


# dictionary with good pairs
global category_per_pair
category_per_pair=symmetrical_hash()
category_per_pair['U']={'A':'a'}
category_per_pair['G']={'A':'k', 'C':'c', 'U':'w'}
category_per_pair_hierarchy='cawko'
# anticodon 1 based position for each model
anticodon = {
    'archaea':38,
    'bacteria':43,
    'eukaryota':38
}

rnaplot_arm_colors = {
    'A':.0,
    'D':.2,
    'C':.4,
    'V':.6,
    'T':.8
    }

def define_extension(hit, args, with_cca = True):
    # define the left and rigth extension of the gene based on the identification of the discriminator base.
    # since the A-T stem in tRNAsec measures 13nt, the discriminator base is always the 14th position starting at base 3' to the T loop
    # if an insertion occurs in the A-T stem, then the 14th nt is already in the infernal hit, in that case the discriminator base is always the base 3' to the prediction

    # check if there is at least one pair in the A stem, otherwise do not extend
    if not [x for x in hit.target_ss_pairs if 'A' in x[-1]]:
        return 0,0
    # make sure if the last pair belongs to the T stem, otherwise do not extend
    if not 'T' in hit.target_ss_pairs[-1][-1]:
        return 0,0

    last_T_stem_pair = hit.target_ss_pairs[-1] # distal pair in the T stem
    i,j = last_T_stem_pair[:2] # given that i and j make a pair, j>i
    AT_length = 13
    if args.AT12:
        AT_length = 12
    # compute the extension necessary to obtain the 14th base from j
    dbase_pos = j + AT_length + 1
    right = max(1, dbase_pos - len(hit.sequence()))+ 3*int(with_cca)
#    right =  max( -(len(hit.sequence()) - j - AT_length) + 1, 1) + 3*int(with_cca)
    left = 0
    if hit.strand == '-':
        left,right = right,left
    return left,right

def hit_nice_output(hit):
    txt_output = '>>tRNA-Sec:%s\n' % str(hit.indexx) 
    txt_output += 'model:%s\n' % hit.query.chromosome
    txt_output += 'chromosome:%s\n' % (hit.chromosome)
    txt_output += 'strand:%s\n' % (hit.strand)
    txt_output += 'positions:%s - %s\n' % (hit.boundaries()[0], hit.boundaries()[1])
    txt_output += 'infernal score:%s\n' % (hit.score)
    txt_output += 'E-value:%s\n' % (hit.evalue)
    txt_output += 'discriminator base:%s\n' % hit.dbase
    txt_output += 'anticodon:%s\n' % hit.acd

    txt_output += '\n\n'
    txt_output += '\t'+hit.target_hand_rf+'\n'
    txt_output += '\t'+hit.extended_seq+'\n'
    txt_output += '\t'+hit.target_ss+'\n'
    txt_output += '\n\n'

    return txt_output

    
def target_ss_pairs(hit):
    target_ss_pairs = []
    for i,j in hit.get_pairs():
        if not hit.hand_rf[i] == hit.hand_rf[j]:
            raise Exception, 'ERROR: wrong hand_rf annotation: %s-%s for this pair: %s-%s' % (hit.hand_rf[i], hit.hand_rf[j], i,j)
        nt_i = hit.alignment.seq_of('t')[i]
        nt_j = hit.alignment.seq_of('t')[j]
        try:
            category=category_per_pair.get(nt_i,nt_j)
        except KeyError:
            category='o'
        i_in_seq = hit.alignment.position_in_seq('t',i+1)-1
        j_in_seq = hit.alignment.position_in_seq('t',j+1)-1
        target_ss_pairs.append( (i_in_seq,j_in_seq,category,hit.hand_rf[i]) )
    return target_ss_pairs

def refine_structure(hit):
    pairs = hit.target_ss_pairs

    # look for good pairs in loops
    pairs_to_be_removed = []
    for n in range( len(hit.target_ss_pairs) ):
        i,j,c,a = hit.target_ss_pairs[n]
        if j-i-1<=3: # a loop shoter than 3
            pairs_to_be_removed.append( hit.target_ss_pairs[n] )
            continue

        try:
            next_pair = hit.target_ss_pairs[n+1]
            next_i = next_pair[0]
        except IndexError:
            # last pair corresponds to the T loop, which also needs to be refined
            next_i = j + 1 # this is to make sure the next if evaluates True

        if j < next_i: # i and j are the loop boundaries
            # check inside the loop
            i_loop = i + 1
            j_loop = j - 1
            while j_loop-i_loop+1>4: # minimum loop of 3 nt
                nt_i = hit.sequence()[i_loop].upper()
                nt_j = hit.sequence()[j_loop].upper()
                try:
                    category=category_per_pair.get(nt_i,nt_j)
                except KeyError:
                    category='o'

                if category in 'ko':
                    break
                pairs.append( (i_loop,j_loop,category,a) )
                hit.refined_ss = hit.ss[:i_loop]+'('+hit.ss[i_loop+1:j_loop]+')'+hit.ss[j_loop+1:]
                i_loop += 1
                j_loop -= 1
    # remove bad pairs identified so far
    for n in pairs_to_be_removed: pairs.remove(n)
    hit.target_ss_pairs = sorted(pairs, key = lambda x: x[0])

    # check if there are unpaired nt between A and T stem, there should not
    a_stem=[]
    t_stem=[]
    for n in hit.target_ss_pairs:
        i,j,c,a = n
        if 'A' in a:
            a_stem.append(n)
        elif 'T' in a:
            t_stem.append(n)
    if not a_stem or not t_stem:
        return
    try:
        last_j_Astem = a_stem[-1][1]
    except:
        print hit
        print a_stem
        print t_stem
        raise

    first_j_Tstem = t_stem[0][1]

    t_stem_length = {'bacteria':5, 'archaea':4, 'eukaryota':4}
    if last_j_Astem - first_j_Tstem > 1:
        # for now it only fixes the case of a single unpaired nt, and a complete T stem
        if last_j_Astem - first_j_Tstem == 2:
            nts_in_t_stem = t_stem[0][1]-t_stem[-1][1]+1
            if len(t_stem) == t_stem_length[ hit.query.chromosome ] or nts_in_t_stem == t_stem_length[ hit.query.chromosome ]: # T stem complete, let's try to extend the A stem
                next_i = a_stem[-1][0]+1
                next_j = a_stem[-1][1]-1
                next_i_nt = hit.sequence()[next_i]
                next_j_nt = hit.sequence()[next_j]
                try:
                    category=category_per_pair.get(nt_i,nt_j)
                except KeyError:
                    category='o'
                if not category in 'ko':
                    pairs.append( (next_i,next_j,category,'A') )
                    hit.refined_ss = hit.ss[:next_i]+'('+hit.ss[next_i+1:next_j]+')'+hit.ss[next_j+1:]
                    hit.target_ss_pairs = sorted(pairs, key = lambda x: x[0])
    


                    
def hit_id(hit):
    if hit.query.chromosome == 'archaea': hit.id = 'tRNA-Sec(a)'
    elif hit.query.chromosome == 'bacteria': hit.id = 'tRNA-Sec(b)'
    elif hit.query.chromosome == 'eukaryota': hit.id = 'tRNA-Sec(e)'


def filter_tRNASec(hit):
    """filter by tRNA-Sec specific features: 
variable arm pairs >= 5 """
    arms_lengths = {}
    for i,j,category,arm in hit.target_ss_pairs:
#        if not category in 'ko': # a -> AT pair    c -> CG pair   k -> GA pair w-> GU pairs   o-> others
        arms_lengths[ arm ] = arms_lengths.setdefault(arm,0) + 1
    # if some arm is missing, skip prediction
    for arm in 'ACDTV':
        if not arm in arms_lengths:
            return False
    # check minimum length of stems
    if arms_lengths['V'] < 5: return False
    if arms_lengths['T'] < 4: return False
    return True


def check_target(infernal_output, target):
    """ """
    for i in open(infernal_output):
        if i.startswith('# target sequence database:'):
            target_infernal = i.rstrip().split()[-1]
            if os.path.abspath( target ) != os.path.abspath( target_infernal ):
                b = bash('md5sum '+ target)
                if b[0]: raise Exception
                target_md5 = b[1].split()[0]
                b = bash('md5sum '+ target_infernal)
                if b[0]: raise Exception, b[1]
                target_infernal_md5 = b[1].split()[0]
                if target_md5 == target_infernal_md5:
                    return True,target_infernal
                else:
                    return False,target_infernal
            elif target == target_infernal:
                return True,target_infernal

def target_ss(hit):
    if not hit.target_ss_pairs: return None
    pairs = {}
    for i,j,c,a in hit.target_ss_pairs:
#        if c in 'ko': continue
        pairs[i]=(j,c,a)
        pairs[j]=(i,c,a)
    ss = ''
    hand_rf = ''
    for n,i in enumerate(hit.sequence()):
        if n in pairs:
            j,c,a = pairs[n]
            if c in 'ko':
                ss+='.'
                hand_rf+='.'
            else:
                if n < pairs[n][0]:
                    ss+='('
                else:
                    ss+=')'
                # hand_rf for this position
                pos_in_ali = hit.alignment.position_in_ali('t',n+1)-1 # 0 based
                hand_rf += hit.hand_rf[pos_in_ali]
 
        else:
            ss+='.'
            pos_in_ali = hit.alignment.position_in_ali('t',n+1)-1 # 0 based
            hand_rf_char = hit.hand_rf[pos_in_ali]
            if hand_rf_char in 'ADCVT':
                hand_rf += '.'
            else:
                hand_rf += hit.hand_rf[pos_in_ali]

    return ss,hand_rf


def pairs_for_pre_string(hit):
    return [(i+1, j+1, category, arm) for i,j,category, arm in hit.target_ss_pairs]


def draw_tRNA(hit, output_folder):

    images_folder = Folder(output_folder+'images/')

    acd = hit.hand_rf.find('acd') + 1 # 1 based
    hit.anticodon_5prime_pos = hit.alignment.position_in_seq('t', acd)
    hit.pairs_for_pre_string = pairs_for_pre_string(hit)

    if len(hit.extended_seq) != len(hit.target_ss):
        raise Exception, 'ERROR: seq lenght:%s and ss length:%s are not the same' % (len(hit.extended_seq), len(hit.target_ss))
    # control paired parentheses
    if hit.target_ss.count('(') != hit.target_ss.count(')'):
        raise Exception, 'ERROR: malformated ss, unbalanced parentheses: >tRNAsec.%s\n%s\n%s' % (hit.indexx,hit.sequence(),hit.target_ss)

    # input for RNAplot
    rnaplot_input =  '>tRNAsec.%s\n%s\n%s' % (hit.indexx, hit.extended_seq, hit.target_ss)

    # generate --pre string for setting the pairs and the colors depending on the arm
    pre = '--pre "/dbasemark { newpath 1 sub coor exch get aload pop fsize 1.75 div 0 360 arc stroke} bind def /outlinecolor {0.6 setgray} bind def /paircolor {0.6 setgray} bind def '
    for i,j,category,arm in hit.pairs_for_pre_string:
        sat=.5
        if not category in 'ko':
            if category ==  'w': sat=.25
            pre +=  ' %s %s %s %s colorpair' % (i,j,rnaplot_arm_colors[arm],sat)
    # anticodon position highlight in pre string
    for i in range(3):
        pre += ' %s cmark' % (hit.anticodon_5prime_pos + i)

    # cmark in discriminator base
    if hit.has_dbase:
        pre += ' %s dbasemark' % (hit.dbase_position)
    pre += '"'

    # run RNAplot
    b = bash('cd '+ images_folder +'; echo "'+ rnaplot_input +'" | RNAplot -t 1 '+ pre +'; cd -')
    if 'ERROR' in b[1]: raise Exception, b[1]
    ps_file = images_folder +'tRNAsec.'+ str(hit.indexx) +'_ss.ps'
    # modify coords for CCA seq
    if hit.has_cca:
        cca_coords(ps_file)
    # expand boundingbox. This solves the problem of having the image cut at the cmarks in the anticodon
    expand_boundingbox(ps_file)
    # convert ps to png
    b = bash('convert -density 150 '+ ps_file +' '+ images_folder +'tRNAsec.'+ str(hit.indexx) +'.png')
    if b[0]: raise Exception, b[1]


def cca_coords(ps_file):
    output = ''
    fh = open(ps_file)
    line = fh.readline()
    output += line
    while line:
        if line.startswith('/sequence'):
            output += line
            line = fh.readline()
            seq = line.rstrip().replace('\\','')
        elif line.startswith('/coor'):
            output += line
            line = fh.readline()
            coords = []
            c = 0
            x = 0
            while not 'def' in line:
                coords.append(tuple(map(float, line.rstrip().strip('[]').split())))
                if c >= len(seq)-3:
                    if not x:
                        x = coords[len(seq)-4][0] + 15
                    else:
                        x += 11
                    y = coords[len(seq)-4][1]
                    line = '[%s %s]\n' % (x, y)
                output += line
                line = fh.readline()
                c += 1
        output += line
        line = fh.readline()
    fh = open(ps_file, 'w')
    print >> fh, output
    fh.close()

def expand_boundingbox(ps_file, expand=10):
    """ """
    output = ''
    for line in open(ps_file):
        if line.startswith('%%BoundingBox:'):
            dsc_header = line.split()[0]
            llx, lly, urx, ury = map(int, line.rstrip().split()[1:])
            lly = lly - expand
            llx = llx - expand
            line = '%s %s %s %s %s\n' % (dsc_header,llx, lly, urx, ury)
        output += line
    fh = open(ps_file, 'w')
    print >> fh, output
    fh.close()

def check_depencies(args):
    # check if cmsearch is available and its version
    b = bash('cmsearch -h')
    if b[0]:
        if 'command not found' in b[1]:
            raise Exception, 'ERROR: cmsearch program not available'
        else:
            raise Exception, b[1]
    for i in b[1].split('\n'):
        if i.startswith('# INFERNAL'):
            if not '# INFERNAL 1.1.1' in i:
                raise Exception, 'ERROR: Infernal version must be 1.1.1, current version is: '+ i
            break
    # check if esl-sfetch 
    b = bash('esl-sfetch -h')
    if b[0]:
        if 'command not found' in b[1]:
            raise Exception, 'ERROR: esl-sfetch program not available'
        else:
            raise Exception, b[1]
    # check if RNAplot is available
    if args.plot:
        b = bash('RNAplot -h')
        if b[0]:
            if 'command not found' in b[1]:
                raise Exception, 'ERROR: RNAplot program not available'
            else:
                raise Exception, b[1]

if __name__=='__main__':
    parser = get_parser()
    args = parser.parse_args()
    main(args)

