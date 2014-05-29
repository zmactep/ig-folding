#!/usr/bin/env python
# -*- coding: utf-8 -*-
# :noTabs=true:

# (c) Copyright Rosetta Commons Member Institutions.
# (c) This file is part of the Rosetta software suite and is made available under license.
# (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
# (c) For more information, see http://www.rosettacommons.org. Questions about this can be
# (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

## @file   antibody.py
## @brief  Pre-processing script for antibody protocol
## @author Sergey Lyskov
## @author Daisuke Kuroda
## @author JKLeman
## @author Jeff Gray

import os
import sys
import re
import json
import subprocess
import shutil

from optparse import OptionParser
from time import time

_script_path_ = os.path.dirname(os.path.realpath(__file__))

_framework_names_ = ['FRL', 'FRH', 'light', 'heavy', 'L1', 'L2', 'L3', 'H1', 'H2', 'H3', 'light_heavy']
_camelid_framework_names_ = ['FRH', 'H1', 'H2', 'H3', 'heavy']

'''
_alignment_legend_to_pretty_legend = {
    'subject-id': 'Subject id',
    'resolution':'Resolution',
    '%-identity': '% identity',
    'alignment-length': 'alignment length',
    'mismatches': 'mismatches',
    'gap-opens': 'gap-opens',
    'q.start': 'q.start',
    'q.end': 'q.end',
    's.start': 's.start',
    's.end': 's.end',
    'evalue': 'evalue',
    'bit-score': 'bit-score'
}'''


def main(args):
    # Script for preparing detecting antibodys and preparing info for Rosetta protocol.
    starttime = time()

    global Filters
    Filters = {filter_by_sequence_length: True, filter_by_alignment_length: False, filter_by_template_resolution: True,
               filter_by_outlier: True, filter_by_template_bfactor: True, filter_by_sequence_homolog: True}

    global Options
    Options = get_options(args, Filters)

    global frlh_info
    global frl_info
    global frh_info

    frlh_info, legend = {}, ''
    for l in open(_script_path_ + '/info/frlh_info'):
        if l.startswith('# '):
            legend = l[2:].split()
        elif len(l) > 8:
            frlh_info[l.split()[0]] = dict(zip(legend, l.split()))

    frl_info, legend = {}, ''
    for l in open(_script_path_ + '/info/frl_info'):
        if l.startswith('# '):
            legend = l[2:].split()
        elif len(l) > 8:
            frl_info[l.split()[0]] = dict(zip(legend, l.split()))

    frh_info, legend = {}, ''
    for l in open(_script_path_ + '/info/frh_info'):
        if l.startswith('# '):
            legend = l[2:].split()
        elif len(l) > 8:
            frh_info[l.split()[0]] = dict(zip(legend, l.split()))

    global script_dir
    script_dir = os.path.dirname(__file__)

    #Options
    Options = process_options(Options)

    if Options.self_test:
        self_test()
        return

    if not (Options.light_chain and Options.heavy_chain or Options.heavy_chain and Options.camelid):
        print('Script for preparing detecting antibodies and preparing info for Rosetta protocol.')
        print('At minimum you need to specify options (--light-chain and --heavy-chain) or (--heavy-chain and --camelid).')
        print('For full list of options run "antibody.py --help"\nERROR: No input chains was specified... exiting...')
        sys.exit(1)

    global framework_names
    if Options.camelid:
        framework_names = _camelid_framework_names_
    else:
        framework_names = _framework_names_

    #read fasta files
    if not Options.camelid:
        light_chain = read_fasta_file(Options.light_chain)
    heavy_chain = read_fasta_file(Options.heavy_chain)

    print('Rosetta Antibody script [Python, version 2.0]. Starting...')

    #create directories
    prefix_details = Options.prefix + 'details/'
    if not os.path.isdir(Options.prefix):
        print('Could not find %s... creating it...' % Options.prefix)
        os.makedirs(Options.prefix)
    if not os.path.isdir(prefix_details):
        print('Could not find %s... creating it...' % prefix_details)
        os.makedirs(prefix_details)

    print('Output prefix:', Options.prefix)
    print('Blast database:', Options.blast_database)
    print('Antibody database:', Options.antibody_database)
    print('rosetta_bin:', Options.rosetta_bin, ' [platform: %s]' % Options.rosetta_platform)
    print('rosetta_database:', Options.rosetta_database)

    #unpacking database files
    antibody_database_files = os.listdir(Options.antibody_database)
    bz2ipped_files = [f for f in antibody_database_files if
                      f.endswith('.bz2') and f[:-4] not in antibody_database_files]
    if bz2ipped_files:
        print('Unpacking rosetta_database files (this need to be done only once)...')
        commandline = 'cd %s && bunzip2 -k %s' % (Options.antibody_database, ' '.join(bz2ipped_files))
        p = subprocess.Popen(commandline, shell=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
        output = p.communicate()[0].decode()
        res = p.poll()

        if Options.verbose and not res:
            print(commandline + '\n', output)
        if res:
            print(commandline + '\n', 'ERROR: Unpacking antibody database files failed with code %s and message: %s' % (
                  res, output))
            sys.exit(1)

    ##### Homolog exclusions settings #####
    global sid_cutoff_cdr
    global sid_cutoff_fr

    if Options.exclude_homologs:  # default exclusion settings if option is set
        if Options.homolog_exclusion == 200:
            Options.homolog_exclusion = 80

    # override default settings
    sid_cutoff = Options.homolog_exclusion
    sid_cutoff_cdr = Options.homolog_exclusion_cdr
    sid_cutoff_fr = Options.homolog_exclusion_fr

    if sid_cutoff_cdr > 100 and sid_cutoff_fr > 100 and sid_cutoff < 100:
        sid_cutoff_cdr = sid_cutoff
        sid_cutoff_fr = sid_cutoff
    if sid_cutoff > 100 and sid_cutoff_cdr > 100 and sid_cutoff_fr > 100:
        print('Using full antibody database (no homolog exclusion)')
    elif sid_cutoff_cdr <= 100 or sid_cutoff_fr <= 100:
        if sid_cutoff_cdr <= 100:
            print('\n!!! Homologs will be excluded with %s SID cut-off during ***CDR*** template selections !!!' \
                  % sid_cutoff_cdr)
        else:
            print('\n!!! Homologs will not be excluded during ***CDR*** template selection (default)    !!!')

        if sid_cutoff_fr <= 100:
            print('\n!!! Homologs will be excluded with %s SID cut-off during ***FR*** template selections !!!' \
                  % sid_cutoff_fr)
        else:
            print('\n!!! Homologs will not be excluded during ***FR*** template selection (default)    !!!')
    ### end Homolog exclusion settings

    if Options.quick:
        Options.relax = 0
        Options.idealize = 0

    print('Idealize:', bool(Options.idealize))
    print('Relax:', bool(Options.relax))

    for name in framework_names:
        if getattr(Options, name):
            print('Custom %s template:' % name, getattr(Options, name))
    print()
    if not Options.camelid:
        print("Light chain: %s" % light_chain)
    print("Heavy chain: %s" % heavy_chain)

    #returns dictionary with CDRs as values
    if not Options.camelid:
        CDRs = IdentifyCDRs(light_chain, heavy_chain)
    else:
        CDRs = Identify_VHH_CDRs(heavy_chain)

    CDRs.update(Extract_FR_CDR_Sequences(**CDRs))
    if not Options.camelid:
        CDRs['light_heavy'] = CDRs['light'] + CDRs['heavy']
    write_results(CDRs, prefix_details)

    if Options.verbose:
        print('CDR:', json.dumps(CDRs, sort_keys=True, indent=2))
    else:
        c = dict(CDRs)
        if not Options.camelid:
            c.pop('numbering_L')
        c.pop('numbering_H')
        # print('CDR:', json.dumps(c, sort_keys=True, indent=2))
    # run Blast
    alignment, legend = run_blast(CDRs, prefix=prefix_details, blast=Options.blast,
                                  blast_database=Options.blast_database, verbose=Options.verbose)

    # create and thread template
    create_virtual_template_pdbs(prefix=prefix_details)
    thread_template_pdbs(CDRs, prefix=prefix_details)

    # superimpose template PDBs
    if not Options.no_superimpose:
        if Options.superimpose_PyRosetta:
            print("\nRunning superimpose_PyRosetta...")
            command = script_dir + "/superimpose_interface.py --prefix " + prefix_details
            p = subprocess.Popen(command, shell=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
            output = p.communicate()[0].decode()
            res = p.poll()

            if Options.verbose and not res:
                print(command + '\n', output)
            if res:
                print(command + '\n', 'ERROR: superimpose_PyRosetta failed.  Code %s\nOutput:\n%s' % (res, output))
                sys.exit(1)
        else:
            print("\nRunning ProFit...")
            superimpose_templates(CDRs, prefix=prefix_details)

    #run Rosetta assemble CDRs
    if Options.rosetta_database:
        run_rosetta(CDRs, prefix=Options.prefix, rosetta_bin=Options.rosetta_bin,
                    rosetta_platform=Options.rosetta_platform, rosetta_database=Options.rosetta_database)
    else:
        print('Rosetta database was not found... skipping rosetta run...')

    results = {'cdr': {}, 'numbering': {}, 'alignment': alignment, 'alignment.legend': legend}
    for k in [i for i in CDRs if not i.startswith('numbering_')]:
        results['cdr'][k] = CDRs[k]
    for n in [i for i in CDRs if i.startswith('numbering_')]:
        results['numbering'][n] = CDRs[n]

    #report run time
    endtime = time()
    print('Run time %2.2f s' % (endtime - starttime))

    with open(Options.prefix + 'results.json', 'w') as f:
        json.dump(results, f, sort_keys=True, indent=2)
    print('Done!')


def get_options(args, filters):
    parser = OptionParser(usage="usage: %prog [OPTIONS] [TESTS]")
    parser.set_description(main.__doc__)

    parser.add_option('-L', '--light-chain',
                      action="store",
                      help="Specify file with the light chain - pure IUPAC ASCII letter sequence, no FASTA headers.",)

    parser.add_option('-H', '--heavy-chain',
                      action="store",
                      help="Specify file with the heavy chain - pure IUPAC ASCII letter sequence, no FASTA headers.",)

    parser.add_option('-c', '--camelid',
                      action="store_true", default=False,
                      help="Handle chain as camelid",)

    parser.add_option('--prefix',
                      action="store", default='grafting/',
                      help="Prefix for output files (directory name). Default is grafting/.",)

    parser.add_option('--blast',
                      action="store", default='blastp',
                      help="Specify path+name for 'blastall' executable. Default is blastp [blast+].",)

    parser.add_option('--superimpose-profit',
                      action="store", default='profit',
                      help="(default).  Add full path as argument if necessary.",)

    parser.add_option('--superimpose-PyRosetta', '--superimpose-pyrosetta',
                      action="store_true", default=False,
                      help="Use PyRosetta to superimpose interface (boolean)",)

    parser.add_option('--no-superimpose', action='store_true', default=False,
                     help='Do not use any superimpose')

    parser.add_option('--blast-database',
                      action="store", default=None,
                      help="Specify path of blast database dir.",)

    parser.add_option('--antibody-database',
                      action="store", default=None,
                      help="Specify path of antibody database dir. Default: script_dir/antibody_database.",)

    parser.add_option('--rosetta-database',
                      action="store", default=None,
                      help="Specify path of rosetta database dir.",)

    parser.add_option('--exclude-homologs', '-x',
                      action="store_true", default=False,
                      help="Exclude homologs with default cutoffs",)

    parser.add_option('--homolog_exclusion',
                      default=200, type="int",
                      help="Specify the cut-off for homolog exclusion during CDR or FR template selections.",)

    parser.add_option('--homolog_exclusion_cdr',
                      default=200, type="int",
                      help="Specify the cut-off for homolog exclusion during CDR template selections.",)

    parser.add_option('--homolog_exclusion_fr',
                      default=200, type="int",
                      help="Specify the cut-off for homolog exclusion during FR template selections.",)

    parser.add_option('--rosetta-bin',
                      action="store", default=None,
                      help="Specify path to 'rosetta/source/bin' dir where antibody_graft', idealize and relax executable expected to be found. Default is '$ROSETTA/main/source/bin, then <script location>/bin' (plasce symlink there) and if not found corresponding steps will be skipped.",)

    parser.add_option('--rosetta-platform',
                      action="store", default=None,
                      help="Specify full extra+compier+build type for rosetta biniaries found in --rosetta-bin. For example use static.linuxgccrelease for static build on Linux. Default is dynamic release build of current OS",)

    parser.add_option("--idealize",
                      action="store_true", default=False, dest="idealize",
                      help="Use idealize protocol on final model.",)

    parser.add_option("--constant-seed",
                      action="store_true", default=False, dest="constant_seed",
                      help="Use constant-seed flag in Rosetta grafting run (for debugging).",)

    parser.add_option("--idealizeoff", "--noidealize",
                      action="store_false", dest="idealize",
                      help="Do not use idealize protocol on final model. (default)",)

    parser.add_option("--relax",
                      action="store_true", default=True, dest="relax",
                      help="Use relax protocol on final model. (default)",)

    parser.add_option("--relaxoff", "--norelax",
                      action="store_false", dest="relax",
                      help="Do not use relax protocol on final model.",)

    parser.add_option("--skip-kink-constraints",
                      action="store_false", dest='kink_constraints', default=True,
                      help="Skip generation of kink constraints file (require PyRosetta). Default is False.",)

    parser.add_option("--timeout",
                      default=900, type="int",
                      help="Maximum runtime for rosetta relax run (use 0 for unlimit), default is 900 - 15min limit",)

    parser.add_option("--quick", "-q",
                      action="store_true", default=False,
                      help="Specify fast run (structure will have clashes).  Prevents stem optimization and turns off "
                           "relax, idealize.",)

    #parser.add_option('--rosetta-options',
    #  action="store", default='',
    #  help="Specify extra options for antibody_graft run.",
    #)

    for name in _framework_names_:
        parser.add_option('--' + name,
                          action="store", default='',
                          help="Specify path or PDB code for %s template. If specified this will overwrite blast "
                               "selection." % name,)

    parser.add_option('--self-test',
                      action="store_true",
                      help="Perform self test by using data in test/ dir and exit.",)

    parser.add_option('--self-test-dir',
                      action="store", default='self-test/',
                      help="Specify path for self test dir [default:self-test/].",)

    parser.add_option('-v', "--verbose",
                      action="store_true", default=False,
                      help="Generate verbose output.",)

    # Filter list of 'filter-function:default state' pairs here, and  extend it to add more filters

    for f in filters:
        parser.add_option('--' + f.__name__.replace('_', '-'), type="int", default=int(filters[f]),
                          help="Boolean option [0/1] that control filtering results with %s function." % f.__name__)

    (options, args) = parser.parse_args(args=args[1:])
    return options

def process_options(options):
    if options.prefix and options.prefix[-1] != '/':
        options.prefix += '/'

    if not options.blast_database:
        options.blast_database = script_dir + '/blast_database'
    if not options.antibody_database:
        options.antibody_database = script_dir + '/antibody_database'

    if not options.rosetta_bin:
        if 'ROSETTA' in os.environ:
            options.rosetta_bin = os.path.abspath(os.environ['ROSETTA']) + '/main/source/bin'
        else:
            options.rosetta_bin = script_dir + '/bin'

    if not options.rosetta_database:
        if os.path.isdir(script_dir + '/rosetta_database'):
            options.rosetta_database = os.path.abspath(script_dir + '/rosetta_database')
        elif 'ROSETTA3_DB' in os.environ:
            options.rosetta_database = os.path.abspath(os.environ['ROSETTA3_DB'])
        elif os.path.isdir(os.environ['HOME'] + '/rosetta_database'):
            options.rosetta_database = os.path.abspath(os.environ['HOME'] + '/rosetta_database')

    options.blast_database = os.path.abspath(options.blast_database)
    options.antibody_database = os.path.abspath(options.antibody_database)
    options.rosetta_bin = os.path.abspath(options.rosetta_bin)

    if not options.rosetta_platform:
        if sys.platform.startswith("linux"):
            options.rosetta_platform = 'linuxgccrelease'
        elif sys.platform == "darwin":
            options.rosetta_platform = 'macosgccrelease'
        else:
            options.rosetta_platform = '_unknown_'

    return options


########################################################
def read_fasta_file(file_name):
    return ''.join([l.rstrip() for l in open(file_name) if not l.startswith('>')]).replace(' ', '').upper()


def write_fasta_file(file_name, data, prefix):
    with open(prefix + file_name + '.fasta', 'w') as f:
        f.write('> %s\n' % file_name)
        f.write(data)
        f.write('\n')


def write_results(CDRs, prefix):
    #with open(prefix+'cdr.json', 'w') as f: json.dump(CDRs, f, sort_keys=True, indent=2)
    for k in [i for i in CDRs if not i.startswith('numbering_')]:
        write_fasta_file(k, CDRs[k], prefix)
    for n in [i for i in CDRs if i.startswith('numbering_')]:
        with open(prefix + n + '.txt', 'w') as f:
            f.write('\n'.join(
                ['%s %s' % (CDRs[n][k], k) for k in sorted(CDRs[n].keys(), key=lambda x: (int_(x), x))]) + '\n')


def safelen(seq):
    return 0 if not seq else len(seq)

class Regions(object):
    PATTERNS = {'VL':{'CDR1':{'pattern':
                              r'C[A-Z]{1,17}(WYL|WLQ|WFQ|WYQ|WYH|WVQ|WVR|WWQ|WVK|WYR|WLL|WFL|WVF|WIQ|WYR|WNQ|WHL|WHQ|WYM|WYY)',
                              'cdr_lower_bound':1,
                              'cdr_upper_bound':-3},
                      'CDR3':{'pattern':r'C[A-Z]{1,15}(L|F|V|S)G[A-Z](G|Y)',
                              'cdr_lower_bound':1,
                              'cdr_upper_bound':-4}},
                'VH':{'CDR1':{'pattern':r'C[A-Z]{1,16}(W)(I|V|F|Y|A|M|L|N|G)(R|K|Q|V|N|C|G)(Q|K|H|E|L|R)',
                              'cdr_lower_bound':4,
                              'cdr_upper_bound':-4},
                      'CDR3':{'pattern':r'C[A-Z]{1,33}(W)(G|A|C)[A-Z](Q|S|G|R)',
                              'cdr_lower_bound':3,
                              'cdr_upper_bound':-4}},
                'VHH':{'CDR1':{'pattern':r'C[A-Z]{1,16}(W)(I|V|F|Y|A|M|L|N|D)(R|K|Q|V|N)(Q|K|H|E|L|R)',
                              'cdr_lower_bound':4,
                              'cdr_upper_bound':-4},
                       'CDR3':{'pattern':r'C[A-Z]{1,27}(W|R)G[A-Z](G|R|D|S)',
                              'cdr_lower_bound':3,
                              'cdr_upper_bound':-4}}}
    SUBCHAIN_BOUNDS = {'VL':{'FIRST_UPPER_BOUND_OF_LONG_CHAIN':65,
                             'SECOND_LOWER_BOUND_OF_LONG_CHAIN':50,
                             'SECOND_SUBCHAIN_LENGTH':75,
                             'SHORT_CHAIN_MAX_LENGTH':120,
                             'FIRST_UPPER_BOUND_OF_SHORT_CHAIN':60,
                             'SECOND_LOWER_BOUND_OF_SHORT_CHAIN':50},
                        'VH':{'FIRST_UPPER_BOUND_OF_LONG_CHAIN':65,
                             'SECOND_LOWER_BOUND_OF_LONG_CHAIN':50,
                             'SECOND_SUBCHAIN_LENGTH':95,
                             'SHORT_CHAIN_MAX_LENGTH':140,
                             'FIRST_UPPER_BOUND_OF_SHORT_CHAIN':60,
                             'SECOND_LOWER_BOUND_OF_SHORT_CHAIN':50},
                       'VHH':{'FIRST_UPPER_BOUND_OF_LONG_CHAIN':65,
                              'SECOND_LOWER_BOUND_OF_LONG_CHAIN':50,
                              'SECOND_SUBCHAIN_LENGTH':95,
                              'SHORT_CHAIN_MAX_LENGTH':140,
                              'FIRST_UPPER_BOUND_OF_SHORT_CHAIN':60,
                              'SECOND_LOWER_BOUND_OF_SHORT_CHAIN':50}}
    FR_CONSTRAINTS = {'VL':{'max_FR1_length':24,
                            'FR4_length':12},
                      'VH':{'max_FR1_length':25,
                            'FR4_length':12},
                      'VHH':{'max_FR1_length':26,
                             'FR4_length':12}}
    CDR_CONSTRAINTS = {'VL':{'CDR2_start_shift':16,
                             'CDR2_length':6},
                       'VH':{'CDR2_start_shift':15,
                             'CDR2_end_shift':33},
                       'VHH':{'CDR2_start_shift':15,
                              'CDR2_end_shift':33}}


    @staticmethod
    def get_subchains(chain, chain_type):
        return (chain[:Regions.SUBCHAIN_BOUNDS[chain_type]['FIRST_UPPER_BOUND_OF_LONG_CHAIN']],
                chain[Regions.SUBCHAIN_BOUNDS[chain_type]['SECOND_LOWER_BOUND_OF_LONG_CHAIN']:
                      Regions.SUBCHAIN_BOUNDS[chain_type]['SECOND_LOWER_BOUND_OF_LONG_CHAIN'] +
                      Regions.SUBCHAIN_BOUNDS[chain_type]['SECOND_SUBCHAIN_LENGTH']]) \
               if len(chain) > Regions.SUBCHAIN_BOUNDS[chain_type]['SHORT_CHAIN_MAX_LENGTH'] else \
               (chain[:Regions.SUBCHAIN_BOUNDS[chain_type]['FIRST_UPPER_BOUND_OF_SHORT_CHAIN']],
                chain[Regions.SUBCHAIN_BOUNDS[chain_type]['SECOND_LOWER_BOUND_OF_SHORT_CHAIN']:])


    @staticmethod
    def get_cdr_by_pattern(chain, chain_type, cdr_name):
        subchain_first, subchain_second = Regions.get_subchains(chain, chain_type)
        res = re.search(Regions.PATTERNS[chain_type][cdr_name]['pattern'], chain)
        return res.group()[Regions.PATTERNS[chain_type][cdr_name]['cdr_lower_bound']:
                           Regions.PATTERNS[chain_type][cdr_name]['cdr_upper_bound']] if res else False


    @staticmethod
    def identify_regions(chain, chain_type):
        FR1=''; FR2=''; FR3=''; FR4=''; CDR1=''; CDR2=''; CDR3=''
        CDR1 = Regions.get_cdr_by_pattern(chain, chain_type, 'CDR1')
        print(chain_type + '_CDR1 detected:', CDR1, '(', safelen(CDR1), 'residues )')
        CDR3 = Regions.get_cdr_by_pattern(chain, chain_type, 'CDR3')
        print(chain_type + '_CDR3 detected:', CDR3, '(', safelen(CDR3), 'residues )')
        if CDR1 and CDR3:
            CDR1_start = chain.index(CDR1)
            CDR1_end = CDR1_start + len(CDR1) - 1
            CDR3_start = chain.index(CDR3)
            CDR3_end = CDR3_start + len(CDR3) - 1
            if CDR1_start > Regions.FR_CONSTRAINTS[chain_type]['max_FR1_length']:
                FR1 = chain[CDR1_start - Regions.FR_CONSTRAINTS[chain_type]['max_FR1_length']:CDR1_start]
            else:
                FR1 = chain[:CDR1_start]
            CDR2_start = CDR1_end + Regions.CDR_CONSTRAINTS[chain_type]['CDR2_start_shift']
            if 'CDR2_length' in Regions.CDR_CONSTRAINTS[chain_type]:
                CDR2_end = CDR2_start + Regions.CDR_CONSTRAINTS[chain_type]['CDR2_length']
            else:
                CDR2_end = CDR3_start - Regions.CDR_CONSTRAINTS[chain_type]['CDR2_end_shift']
            CDR2 = chain[CDR2_start: CDR2_end + 1]
            print(chain_type + '_CDR2 detected:', CDR2, '(', safelen(CDR2), 'residues )')

            FR2 = chain[CDR1_end + 1: CDR2_start]
            FR3 = chain[CDR2_end + 1: CDR3_start]
            FR4 = chain[CDR3_end + 1: CDR3_end + Regions.FR_CONSTRAINTS[chain_type]['FR4_length'] + 1]

            print(chain_type + '_FR1:', FR1)
            print(chain_type + '_FR2:', FR2)
            print(chain_type + '_FR3:', FR3)
            print(chain_type + '_FR4:', FR4)
            print(chain_type,'segments:', FR1, CDR1, FR2, CDR2, FR3, CDR3, FR4)
        return {'FR1':FR1, 'FR2':FR2, 'FR3':FR3, 'FR4':FR4, 'CDR1':CDR1, 'CDR2':CDR2, 'CDR3':CDR3}


def IdentifyCDRs(light_chain, heavy_chain):
    L = Regions.identify_regions(light_chain, 'VL')
    H = Regions.identify_regions(heavy_chain, 'VH')

    if not (L['CDR1'] and L['CDR3'] and H['CDR1'] and H['CDR3']):
        if not L['CDR1']:
            print('ERROR: CDR L1 cannot be recognized !!!  L1 pattern: ', Regions.PATTERNS['VL']['CDR1']['pattern'])

        if not L['CDR3']:
            print('ERROR: CDR L3 cannot be recognized !!!  L3 pattern: ', Regions.PATTERNS['VL']['CDR3']['pattern'])

        if not H['CDR1']:
            print('ERROR: CDR H1 cannot be recognized !!!  H1 pattern: ', Regions.PATTERNS['VH']['CDR1']['pattern'])

        if not H['CDR3']:
            print('ERROR: CDR H3 cannot be recognized !!!  H3 pattern: ', Regions.PATTERNS['VH']['CDR3']['pattern'])
        sys.exit(1)

    return(dict(L1=L['CDR1'],L2=L['CDR2'],L3=L['CDR3'],FR_L1=L['FR1'],FR_L2=L['FR2'],FR_L3=L['FR3'],FR_L4=L['FR4'],
                H1=H['CDR1'],H2=H['CDR2'],H3=H['CDR3'],FR_H1=H['FR1'],FR_H2=H['FR2'],FR_H3=H['FR3'],FR_H4=H['FR4']))


def Identify_VHH_CDRs(heavy_chain):
    H = Regions.identify_regions(heavy_chain, 'VHH')

    if not (H['CDR1'] and H['CDR3']):
        if not H['CDR1']:
            print('ERROR: CDR H1 cannot be recognized !!!  H1 pattern: ', Regions.PATTERNS['VHH']['CDR1']['pattern'])

        if not H['CDR3']:
            print('ERROR: CDR H3 cannot be recognized !!!  H3 pattern: ', Regions.PATTERNS['VHH']['CDR3']['pattern'])
        sys.exit(1)

    return dict(H1=H['CDR1'],H2=H['CDR2'],H3=H['CDR3'],FR_H1=H['FR1'],FR_H2=H['FR2'],FR_H3=H['FR3'],FR_H4=H['FR4'])

def int_(s):
    return int(re.sub('[A-Z]', '', s))
    # v = int(re.sub('[A-Z]', '', new_number_FR_L1) )
    #   $new_number_FR_L1[$i] =~ s/[A-Z]//
    # new_number_FR_L1[i] = string.translate(new_number_FR_L1[i], None, string.ascii_letters)
    
class Cut_sequence(object):
    CUT_BEFORE = {'VL':{'FR1_WITH_G':{'max_length':23,'max_cut':1},
                        'FR1_MISSING_G':{'max_length':22,'max_cut':2}},
                  'VH':{'FR1':{'max_length':25,'max_cut':1}}}

    @staticmethod
    def cut_before_if_needed(region, chain_type, region_name):
        d = len(region) - Cut_sequence.CUT_BEFORE[chain_type][region_name]['max_length']
        if d > Cut_sequence.CUT_BEFORE[chain_type][region_name]['max_cut']:
            print('%s of %s has wrong length' % region_name, chain_type)
            sys.exit(1)
        if d > 0:
            return region[d:]
        return region


class Numbering(object):
    REGION_LENGTH_BOUNDS = {'VL':{'FR1':[19,24],
                                  'FR2':[15,15],
                                  'FR3':[32,34],
                                  'FR4':[0,1000],
                                  'CDR1':[8,17],
                                  'CDR2':[7,11],
                                  'CDR3':[5,15]},
                            'VH':{'FR1':[16,26],
                                  'FR2':[14,14],
                                  'FR3':[30,32],
                                  'FR4':[0,1000],
                                  'CDR1':[6,13],
                                  'CDR2':[12,22],
                                  'CDR3':[3,34]}}

    SPECIAL_PATTERNS = {'VL':{'FR1_WITH_G':r'[A-Z][QE][A-Z]{9}[A-Z][A-Z]{4}[LVIMF][A-Z]C',
                              'FR1_MISSING_G':r'[A-Z][QE][A-Z]{8}[A-Z][A-Z]{4}[LVIMF][A-Z]C'}}

    # numbering pattern:
    # ranges of numbers which are used for  numbering are divided by ','.
    # bounds of ranges are splitted by ':'
    # 'A:B,C:D,...' - [A,B] - first range, [C,D] - second range, etc
    # some ranges can be replaced by single number such as 's82A', where s - special symbol and 82A - number for that position
    # One bound of one range can be writen as 'a1-a2*a3'
    # '-' sign mean that bound can expand in band [a1,a2] to match the pattern
    # '*' sign mean which number a3 will be suplemented with letter during expansion if bound range will be overflowed
    NUMBERING_SCHEME = {'VL':{'FR1_WITH_G':'1-5*0:23',
                              'FR1_MISSING_G':'1-5*0:8,10:23',
                              'FR2':'35:49',
                              'FR3':'57:65,66-66*66:88',
                              'FR4':'98:109',
                              'CDR1':'24:30,31-34*30:34',
                              'CDR2':'50:54,55-56*54:56',
                              'CDR3':'89:92-96*95,97:97'},
                        'VH':{'FR1':'1-10*0:25',
                              'FR2':'36:49',
                              'FR3':'66:73-75*75,76:82,s82A,s82B,s82C,83:94',
                              'FR4':'103:114',
                              'CDR1':'26:27-31*31,32:35',
                              'CDR2':'50:52,53-57*52:65',
                              'CDR3':'95:97-102*100'}}

    class EXPANSION_FLAG(object):
        LEFT = 1
        RIGHT = 2


    @staticmethod
    def get_numbers(l,r):
        return ','.join(str(x) for x in range(l,r+1))


    @staticmethod
    def get_free_bound_numbers(b_fixed, b_min, b_max, b_ext, diff, expansion_flag):
        extra_numbering = ''
        if expansion_flag == Numbering.EXPANSION_FLAG.LEFT and b_fixed - b_min + 1 >= diff \
        or expansion_flag == Numbering.EXPANSION_FLAG.RIGHT and b_max - b_fixed + 1 >= diff:
            if expansion_flag == Numbering.EXPANSION_FLAG.LEFT:
                l = b_fixed - diff + 1
                r = b_fixed
            elif expansion_flag == Numbering.EXPANSION_FLAG.RIGHT:
                l = b_fixed
                r = b_fixed + diff - 1
            return Numbering.get_numbers(l,r)
        else:
            if expansion_flag == Numbering.EXPANSION_FLAG.LEFT:
                letters_count = diff - (b_fixed - b_min + 1)
                l = b_min
                r = b_fixed
            elif expansion_flag == Numbering.EXPANSION_FLAG.RIGHT:
                letters_count = diff - (b_max - b_fixed + 1)
                l = b_fixed
                r = b_max
            extra_numbering = ','.join(str(b_ext) + chr(x) for x in range(ord('A'), ord('A') + letters_count))
        return ','.join(filter(None,[Numbering.get_numbers(l,b_ext), extra_numbering, Numbering.get_numbers(b_ext+1,r)]))


    @staticmethod
    def enumerate_chain(region_chain, chain_type, region_name):
        numbering_scheme = Numbering.NUMBERING_SCHEME[chain_type][region_name]
        subranges = numbering_scheme.split(',')
        numberings = [''] * len(subranges)
        symbols_processed = 0
        subrange_with_free_bound = -1
        extra_numbering = ''
        for i in range(len(subranges)):
            if subranges[i][0] == 's':
                numberings[i] = subranges[i][1:]
                symbols_processed += 1
            elif '-' in subranges[i]:
                subrange_with_free_bound = i
            else:
                l,r = [int(x) for x in subranges[i].split(':')]
                numberings[i] = Numbering.get_numbers(l,r)
                symbols_processed += r - l + 1
        if subrange_with_free_bound != -1:
            l,r = subranges[subrange_with_free_bound].split(':')
            diff = len(region_chain) - symbols_processed
            if '-' in l:
                band = l
                b_fixed = int(r)
                expansion_flag = Numbering.EXPANSION_FLAG.LEFT
            elif '-' in r:
                b_fixed = int(l)
                band = r
                expansion_flag = Numbering.EXPANSION_FLAG.RIGHT
            b_min, b_max, b_ext = [int(x) for x in sum([x.split('*') for x in band.split('-')],[])]
            numberings[subrange_with_free_bound] = Numbering.get_free_bound_numbers(b_fixed, b_min, b_max, b_ext, diff, expansion_flag)
        return ','.join(filter(None,numberings))


    @staticmethod
    def check_region_length(region_sequence, chain_name, region_name):
        return (len(region_sequence) >= Numbering.REGION_LENGTH_BOUNDS[chain_name][region_name][0] and
                len(region_sequence) <= Numbering.REGION_LENGTH_BOUNDS[chain_name][region_name][1])


    @staticmethod
    def process_region(region_sequence, chain_name, region_name):
        if not Numbering.check_region_length(region_sequence, chain_name, region_name):
            print('ERROR: wrong length of ' + region_name + ' region in ' + chain_name + ' sequence and is not between [' +
                  str(Numbering.REGION_LENGTH_BOUNDS[chain_name][region_name][0]) + ','
                  + str(Numbering.REGION_LENGTH_BOUNDS[chain_name][region_name][1]) + ']')
            exit(1)
        if region_name == 'FR1':
            if chain_name == 'VL':
                if re.search(Numbering.SPECIAL_PATTERNS['VL']['FR1_WITH_G'], region_sequence):
                    region_name = 'FR1_WITH_G'
                elif re.search(Numbering.SPECIAL_PATTERNS['VL']['FR1_MISSING_G'], region_sequence):
                    region_name = 'FR1_MISSING_G'
                else:
                    print('ERROR: Current code could not assign Chothia numbering of VL:FR1 in the query sequence!!! Exiting...')
                    exit(1)
            region_sequence = Cut_sequence.cut_before_if_needed(region_sequence,chain_name,region_name)
        new_number = Numbering.enumerate_chain(region_sequence,chain_name,region_name)

        return region_sequence, new_number



def Extract_FR_CDR_Sequences(L1='', L2='', L3='', H1='', H2='', H3='', FR_L1='', FR_L2='', FR_L3='', FR_L4='', FR_H1='',
                             FR_H2='', FR_H3='', FR_H4=''):
    camelid = False
    if Options.camelid:
        camelid = True
    
    if not camelid:
        regions = {'VL':{'CDR1':L1,'CDR2':L2,'CDR3':L3,'FR1':FR_L1,'FR2':FR_L2,'FR3':FR_L3,'FR4':FR_L4},
                   'VH':{'CDR1':H1,'CDR2':H2,'CDR3':H3,'FR1':FR_H1,'FR2':FR_H2,'FR3':FR_H3,'FR4':FR_H4}}
        new_number = {'VL':{},'VH':{}}
    else:
        regions = {'VH':{'CDR1':H1,'CDR2':H2,'CDR3':H3,'FR1':FR_H1,'FR2':FR_H2,'FR3':FR_H3,'FR4':FR_H4}}
        new_number = {'VH':{}}
    
    for chain_name in regions.keys():
        for region_name in regions[chain_name].keys():
            regions[chain_name][region_name], \
            new_number[chain_name][region_name] = Numbering.process_region(regions[chain_name][region_name],
                    chain_name, region_name)

    if not camelid:
        new_number_L1 = new_number['VL']['CDR1']
        new_number_L2 = new_number['VL']['CDR2']
        new_number_L3 = new_number['VL']['CDR3']
        new_number_FR_L1 = new_number['VL']['FR1']
        new_number_FR_L2 = new_number['VL']['FR2']
        new_number_FR_L3 = new_number['VL']['FR3']
        new_number_FR_L4 = new_number['VL']['FR4']

    new_number_H1 = new_number['VH']['CDR1']
    new_number_H2 = new_number['VH']['CDR2']
    new_number_H3 = new_number['VH']['CDR3']
    new_number_FR_H1 = new_number['VH']['FR1']
    new_number_FR_H2 = new_number['VH']['FR2']
    new_number_FR_H3 = new_number['VH']['FR3']
    new_number_FR_H4 = new_number['VH']['FR4']

    if not camelid:
        L1 = regions['VL']['CDR1']
        L2 = regions['VL']['CDR2']
        L3 = regions['VL']['CDR3']
        FR_L1 = regions['VL']['FR1']
        FR_L2 = regions['VL']['FR2']
        FR_L3 = regions['VL']['FR3']
        FR_L4 = regions['VL']['FR4']

    H1 = regions['VH']['CDR1']
    H2 = regions['VH']['CDR2']
    H3 = regions['VH']['CDR3']
    FR_H1 = regions['VH']['FR1']
    FR_H2 = regions['VH']['FR2']
    FR_H3 = regions['VH']['FR3']
    FR_H4 = regions['VH']['FR4']
        
    # Converting all new_number_* vars in to a lists
    if not camelid:
        new_number_L1 = new_number_L1.split(',')
        new_number_L2 = new_number_L2.split(',')
        new_number_L3 = new_number_L3.split(',')
        new_number_FR_L1 = new_number_FR_L1.split(',')
        new_number_FR_L2 = new_number_FR_L2.split(',')
        new_number_FR_L3 = new_number_FR_L3.split(',')
        new_number_FR_L4 = new_number_FR_L4.split(',')
    new_number_H1 = new_number_H1.split(',')
    new_number_H2 = new_number_H2.split(',')
    new_number_H3 = new_number_H3.split(',')
    new_number_FR_H1 = new_number_FR_H1.split(',')
    new_number_FR_H2 = new_number_FR_H2.split(',')
    new_number_FR_H3 = new_number_FR_H3.split(',')
    new_number_FR_H4 = new_number_FR_H4.split(',')

    print_seq_FR_L1, print_seq_FR_L2, print_seq_FR_L3, print_seq_FR_L4, print_seq_FR_L4_extra = '', '', '', '', ''
    print_seq_FR_H1, print_seq_FR_H2, print_seq_FR_H3, print_seq_FR_H4, print_seq_FR_H4_extra = '', '', '', '', ''
    numbering_L, numbering_H = {}, {}

    if not camelid:
    # OUTPUT FOR LIGHT CHAIN. This should be save in the file 'numbering_L.txt'.
        for i, s in enumerate(FR_L1):  # FR_L1
            numbering_L[new_number_FR_L1[i]] = s  # +='%s %s\n' % (s, new_number_FR_L1[i])
            v = int_(new_number_FR_L1[i])
            if 10 <= v <= 23:
                print_seq_FR_L1 += s

        for i, s in enumerate(L1): numbering_L[new_number_L1[i]] = s  # +='%s %s\n' % (s, new_number_L1[i])  # L1

        for i, s in enumerate(FR_L2):  # FR_L2
            numbering_L[new_number_FR_L2[i]] = s  # +='%s %s\n' % (s, new_number_FR_L2[i])
            v = int_(new_number_FR_L2[i])
            if (35 <= v <= 38) or (45 <= v <= 49):
                print_seq_FR_L2 += s

        for i, s in enumerate(L2):
            numbering_L[new_number_L2[i]] = s  # +='%s %s\n' % (s, new_number_L2[i])  # L2

        for i, s in enumerate(FR_L3):  # FR_L3
            numbering_L[new_number_FR_L3[i]] = s  # +='%s %s\n' % (s, new_number_FR_L3[i])
            v = int_(new_number_FR_L3[i])
            if (57 <= v <= 66) or (71 <= v <= 88):
                print_seq_FR_L3 += s

        for i, s in enumerate(L3): numbering_L[new_number_L3[i]] = s  # +='%s %s\n' % (s, new_number_L3[i])  # L3

        for i, s in enumerate(FR_L4):  # FR_L4
            numbering_L[new_number_FR_L4[i]] = s  # +='%s %s\n' % (s, new_number_FR_L4[i])
            v = int_(new_number_FR_L4[i])
            if 98 <= v <= 104:
                print_seq_FR_L4 += s
            if 98 <= v <= 101:
                print_seq_FR_L4_extra += s

    #  OUTPUT FOR HEAVY CHAIN. This should be save in the file 'numbering_H.txt'.
    for i, s in enumerate(FR_H1):  # FR_H1
        numbering_H[new_number_FR_H1[i]] = s  # +='%s %s\n' % (s, new_number_FR_H1[i])
        v = int_(new_number_FR_H1[i])
        if 10 <= v <= 25:
            print_seq_FR_H1 += s

    for i, s in enumerate(H1):
        numbering_H[new_number_H1[i]] = s  #+='%s %s\n' % (s, new_number_H1[i])  # H1

    for i, s in enumerate(FR_H2):  # FR_H2
        numbering_H[new_number_FR_H2[i]] = s  # +='%s %s\n' % (s, new_number_FR_H2[i])
        v = int_(new_number_FR_H2[i])
        if (36 <= v <= 39) or (46 <= v <= 49):
            print_seq_FR_H2 += s

    for i, s in enumerate(H2):
        numbering_H[new_number_H2[i]] = s  # +='%s %s\n' % (s, new_number_H2[i])  # H2

    for i, s in enumerate(FR_H3):  # FR_H3
        numbering_H[new_number_FR_H3[i]] = s  # +='%s %s\n' % (s, new_number_FR_H3[i])
        v = int_(new_number_FR_H3[i])
        if 66 <= v <= 94:
            print_seq_FR_H3 += s

    for i, s in enumerate(H3):
        numbering_H[new_number_H3[i]] = s  # +='%s %s\n' % (s, new_number_H3[i])  # H3

    for i, s in enumerate(FR_H4):  # FR_H4
        numbering_H[new_number_FR_H4[i]] = s  # +='%s %s\n' % (s, new_number_FR_H4[i])
        v = int_(new_number_FR_H4[i])
        if 103 <= v <= 109:
            print_seq_FR_H4 += s
        if 103 <= v <= 106:
            print_seq_FR_H4_extra += s
    if not camelid:
        FRL = print_seq_FR_L1 + print_seq_FR_L2 + print_seq_FR_L3 + print_seq_FR_L4
    FRH = print_seq_FR_H1 + print_seq_FR_H2 + print_seq_FR_H3 + print_seq_FR_H4

    if not camelid and len(FRL) != 58 and len(FRL) != 60:
        print("ERROR: Current DB does not cover the length of FRL of your query.")
        print("ERROR: FRL length of your query:", len(FRL))
        print("ERROR: DB: 58 or 60")
        sys.exit(1)

    if len(FRH) != 63 and len(FRH) != 65:
        print("ERROR: Current DB does not cover the length of FRL of your query.")
        print("ERROR: FRH length of your query:", len(FRH))
        print("ERROR: DB: 63 or 65")
        sys.exit(1)

    if not camelid:
        result = dict(FRL=print_seq_FR_L1 + print_seq_FR_L2 + print_seq_FR_L3 + print_seq_FR_L4,
                FRH=print_seq_FR_H1 + print_seq_FR_H2 + print_seq_FR_H3 + print_seq_FR_H4,
                light=FR_L1 + L1 + FR_L2 + L2 + FR_L3 + L3 + print_seq_FR_L4_extra,
                heavy=FR_H1 + H1 + FR_H2 + H2 + FR_H3 + H3 + print_seq_FR_H4_extra,
                numbering_L=numbering_L, numbering_H=numbering_H)
    else:
        result = dict(FRH=print_seq_FR_H1 + print_seq_FR_H2 + print_seq_FR_H3 + print_seq_FR_H4,
                heavy=FR_H1 + H1 + FR_H2 + H2 + FR_H3 + H3 + print_seq_FR_H4_extra,
                numbering_H=numbering_H)
    return result


def run_blast(cdr_query, prefix, blast, blast_database, verbose=False):
    camelid = False
    if Options.camelid:
        camelid = True
    print('\nRunning %s' % blast)
    cdr_info, legend = {}, ''  # first reading cdr_info table
    for l in open(_script_path_ + '/info/cdr_info'):
        if l.startswith('# '):
            legend = l[2:].split()
        elif len(l) > 8:
            cdr_info[l.split()[0]] = dict(zip(legend, l.split()))

    for line in open(_script_path_ + '/info/fv_length'):
        cdr_info[line.split()[0]][line.split()[1] + '_length'] = int(line.split()[2])

    # cdr_info consistency check
    if not camelid:
        cdrs_names = ['L1', 'L2', 'L3', 'H1', 'H2', 'H3']
    else:
        cdrs_names = ['H1', 'H2', 'H3']
    for name, i in cdr_info.items():
        if sum([i[k] != 'none' and int(i[k + '_length']) != len(i[k]) for k in cdrs_names]):
            print('ERROR: cdr_info length info is inconsistent for line: %s' % name)
            for k in cdrs_names:
                print(k, i[k], i[k + '_length'], len(i[k]))
            sys.exit(1)

    alignment = {'summary': []}
    for k in framework_names:
        input_query = k + '.fasta'
        output = k + '.align'  # Options.prefix +

        # if CDRs, use length-depend DBs for BLAST
        len_cdr = len(cdr_query[k])
        if k.count('FR') or k.count('heavy') or not camelid and k.count('light'):
            db = blast_database + '/database.%s' % k
        else:
            db = blast_database + '/database.%s.%s' % (k, len_cdr)

        # check that database file exists
        if not os.path.isfile(db):
            print('\nERROR: No %s templates of length %s (%s)\n' % (k, len_cdr, db))
            sys.exit(1)

        wordsize, e_value, matrix = (2, 0.00001, 'BLOSUM62') if k.count('FR') or k.count('heavy')\
        or not camelid and k.count('light') else (2, 2000, 'PAM30')

        commandline = 'cd %s && ' % '\ '.join(os.path.dirname(prefix).split()) + blast + \
                      (' -db %(db)s -query %(input_query)s -out %(output)s -evalue %(e_value)s -matrix %(matrix)s' +
                       ' -word_size %(wordsize)s -outfmt 7 -max_target_seqs 600') % vars()

        p = subprocess.Popen(commandline, shell=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
        output = p.communicate()[0].decode()
        res = p.poll()

        if verbose and not res:
            print(commandline + '\n', output)
        if res:
            print(commandline + '\n', 'ERROR: Blast execution failed with code %s and message: %s' % (res, output))
            sys.exit(1)

        # print('Filtering results...')
        table, legend = [], ''
        for l in open(prefix + k + '.align'):
            if l.startswith('# Fields: '):
                legend = [i.strip().replace('. ', '.').replace(' ', '-') for i in l[10:-1].split(',')]
            elif l.startswith('Query_1') or l.startswith(k):
                table.append(dict(zip(legend, l.split())))
                # Sergey: Blast+ 2.2.27 have lines starts with Query_1 and 2.2.28 with name of the file (framework-name)

        table_original = table[:]

        def sort_and_write_results(table, file_name):
            table.sort(key=lambda x: (-float(x['bit-score']), float(
                cdr_info[x['subject-id']]['resolution'])))  # sort results first by bit-score then by resolution
            for r in table:
                r['resolution'] = cdr_info[r['subject-id']]['resolution']
            with open(file_name, 'w') as f:
                f.write('# Fields: ' + ' '.join(legend) + '\n')
                for i in table:
                    f.write('\t'.join([i[k] for k in legend]) + '\n')

                def try_to_float(v):
                    try:
                        return int(v)
                    except ValueError:
                        try:
                            return float(v)
                        except ValueError:
                            return v

                for i, v in enumerate(table):
                    table[i] = dict([(k, try_to_float(v[k])) for k in v])

        prefix_filters = prefix + 'filters/'
        if not os.path.isdir(prefix_filters):
            print('Could not find %s... creating it...' % prefix_filters)
            os.makedirs(prefix_filters)

        for f in Filters:
            if getattr(Options, f.__name__):
                f(k, table, cdr_query, cdr_info)
                t = table_original[:]
                f(k, t, cdr_query, cdr_info)
                sort_and_write_results([i for i in table_original if i not in t],
                                       prefix_filters + 'filtered-by.' + f.__name__[7:] + '.' + k + '.align')

        sort_and_write_results(table, prefix_filters + 'filtered.' + k + '.align')  # Writing filtered results.
        alignment[k] = table
        if table:
            alignment['summary'].append(dict(table[0]))
            alignment['summary'][-1]['subject-id'] = k

        custom_template = getattr(Options, k)

        if custom_template and not os.path.isfile(custom_template):
            custom_template = '/pdb%s_chothia.pdb' % custom_template
        if not custom_template:
            if table:
                custom_template = table[0]['subject-id']
            else:  # if there is no template... table is a list, which has a blast result
                for v in cdr_info.items():
                    check_length = '%s_length' % k
                    if len_cdr == int(v[1][check_length]):
                        pdb_random = v[0]
                        break
                print('\nWARNING: No template avaliable for %s after filtering! Using a random template of the same ' \
                      'length as the query\n' % k)
                custom_template = pdb_random
                # sys.exit(1)
            print("%s template: %s" % (k, custom_template))
        else:
            print('Custom %s template: %s...' % (k, custom_template))
        shutil.copy(Options.antibody_database + '/' + custom_template, prefix + '/template.' + k + '.pdb')

    legend.remove('query-id')
    legend.insert(1, 'resolution')
    return alignment, legend


def create_virtual_template_pdbs(prefix):
    camelid = False
    if Options.camelid:
        camelid = True
    # Make a template PDB file, which has psuedo atoms, so that the query and a template can have the same length sequence.
    if camelid:
        chains_names = ['H']
    else:
        chains_names = ['L', 'H']
    for chain in chains_names:
        with open(prefix + '/template.tmp.FR' + chain + '.pdb', 'w') as o:
            # Count the number of lines of a template PDB file
            cnt1 = 0
            for line2 in open(prefix + '/template.FR' + chain + '.pdb'):
                if line2[0:4] == 'ATOM':
                    chain_temp = line2[21:22]
                    if chain_temp == chain:
                        cnt1 += 1

            cnt_res = 0
            for line in open(prefix + '/numbering_' + chain + '.txt'):
                cnt_res += 1
                check = 0
                res_num_q = line[2:].rstrip('\n')
                cnt2 = 0
                for line2 in open(prefix + '/template.FR' + chain + '.pdb'):
                    if line2[0:4] == 'ATOM':
                        chain_temp = line2[21:22]
                        res_aa = line2[17:20]
                        if chain_temp == chain:
                            cnt2 += 1
                        res_num_temp = line2[23:27].strip()
                        iCode_tmp = line2[27:27]
                        x_coord = line2[31:38]
                        y_coord = line2[39:46]
                        z_coord = line2[47:54]

                        # Identify the residue number of the last residue
                        if cnt2 == 1 and chain_temp == chain:
                            res_num_temp_first = res_num_temp
                        elif cnt2 == cnt1 and chain_temp == chain:
                            res_num_temp_last = res_num_temp

                        if chain == chain_temp and res_num_q == res_num_temp:
                            check = 1
                            o.write(line2)

                # v needed to circumvent collinearity!
                x_coord_n = '%8.3f' % (float(x_coord) + 100 + cnt_res * 2)
                x_coord_ca = '%8.3f' % (float(x_coord) + 110 + cnt_res * 2)
                x_coord_c = '%8.3f' % (float(x_coord) + 120 + cnt_res * 2)
                x_coord_o = '%8.3f' % (float(x_coord) + 130 + cnt_res * 2)

                y_coord_n = '%8.3f' % (float(y_coord) + 100 + cnt_res * 2)
                y_coord_ca = '%8.3f' % (float(y_coord) + 110 + cnt_res * 2)
                y_coord_c = '%8.3f' % (float(y_coord) + 120 + cnt_res * 2)
                y_coord_o = '%8.3f' % (float(y_coord) + 130 + cnt_res * 2)

                z_coord_n = '%8.3f' % (float(z_coord) + 100 + cnt_res * 2)
                z_coord_ca = '%8.3f' % (float(z_coord) + 110 + cnt_res * 2)
                z_coord_c = '%8.3f' % (float(z_coord) + 120 + cnt_res * 2)
                z_coord_o = '%8.3f' % (float(z_coord) + 130 + cnt_res * 2)

                # x_coord_n  = '%8.3f' % (float(x_coord) + random.uniform(1, 200))
                # ...
                # z_coord_o  = '%8.3f' % (float(z_coord) + random.uniform(1, 200))

                if res_num_q.isdigit():
                    if int(res_num_q) < 100:
                        res_num_q_tmp = '  ' + res_num_q
                    elif int(res_num_q) >= 100:
                        res_num_q_tmp = ' ' + res_num_q
                    iCode_q_tmp = ' '
                else:
                    res_num_q_tmp = res_num_q[:-1]
                    if int(res_num_q_tmp) < 100:
                        res_num_q_tmp = '  ' + res_num_q_tmp
                    elif int(res_num_q_tmp) >= 100:
                        res_num_q_tmp = ' ' + res_num_q_tmp
                    iCode_q_tmp = res_num_q[-1:]

                if check == 0 and (int(res_num_temp_first) <= int(res_num_q_tmp) <= int(res_num_temp_last)):
                    r_n = dict(tempFactor=' 25.00', chainID=chain, name=' N  ', altLoc=' ', occupancy='  1.00',
                               element='  ', resSeq=res_num_q_tmp, charge='  ', iCode=iCode_q_tmp, resName='GLY',
                               x=x_coord_n, serial=' 0000', z=z_coord_n, type='ATOM  ', y=y_coord_n)
                    r_ca = dict(tempFactor=' 25.00', chainID=chain, name=' CA ', altLoc=' ', occupancy='  1.00',
                                element='  ', resSeq=res_num_q_tmp, charge='  ', iCode=iCode_q_tmp, resName='GLY',
                                x=x_coord_ca, serial=' 0000', z=z_coord_ca, type='ATOM  ', y=y_coord_ca)
                    r_c = dict(tempFactor=' 25.00', chainID=chain, name=' C  ', altLoc=' ', occupancy='  1.00',
                               element='  ', resSeq=res_num_q_tmp, charge='  ', iCode=iCode_q_tmp, resName='GLY',
                               x=x_coord_c, serial=' 0000', z=z_coord_c, type='ATOM  ', y=y_coord_c)
                    r_o = dict(tempFactor=' 25.00', chainID=chain, name=' O  ', altLoc=' ', occupancy='  1.00',
                               element='  ', resSeq=res_num_q_tmp, charge='  ', iCode=iCode_q_tmp, resName='GLY',
                               x=x_coord_o, serial=' 0000', z=z_coord_o, type='ATOM  ', y=y_coord_o)

                    o.write(records_to_pdb_string(r_n) + '\n')
                    o.write(records_to_pdb_string(r_ca) + '\n')
                    o.write(records_to_pdb_string(r_c) + '\n')
                    o.write(records_to_pdb_string(r_o) + '\n')


def thread_template_pdbs(CDRs, prefix):
    camelid = False
    if Options.camelid:
        camelid = True

    L_regions = dict(L1=(24, 34), L2=(50, 56), L3=(89, 97))
    H_regions = dict(H1=(26, 35), H2=(50, 65), H3=(95, 102))

    if not camelid:
        region_info = [('L', 'FRL', 'numbering_L', L_regions), ('H', 'FRH', 'numbering_H', H_regions)]
    else:
        region_info = [('H', 'FRH', 'numbering_H', H_regions)]

    for chain, k, numbering, regions in region_info:
        with open(prefix + '/template.threaded.' + k + '.pdb', 'w') as o:
            for line in open(prefix + '/template.tmp.' + k + '.pdb'):  #  This is a template PDB in our database
                #for line in open(prefix+'/template.'+k+'.pdb'): # This is a template PDB in our database
                r = map_pdb_string_to_records(line)  # 'r' is a info of a template
                if r['type'] == "ATOM  ":
                    res_num = int(r['resSeq'])
                    res_num_icode = ('%s%s' % (res_num, r['iCode'])).strip()
                    if res_num_icode in CDRs[numbering] and chain == r['chainID']:
                        for loop in regions:
                            if regions[loop][0] <= res_num <= regions[loop][1]:
                                # r['charge'] = '.' if r['resName'] == AA_Code[ CDRs[numbering][res_num_icode] ] else '!'
                                if r['name'] in [' N  ', ' CA ', ' C  ', ' O  ']:
                                # r['resName'] == AA_Code[ CDRs[numbering][res_num_icode] ]
                                    for coord in 'xyz':
                                        r[coord] = '%8.3f' % (float(r[coord]) + 100.0)
                                    r['resName'] = AA_Code[CDRs[numbering][res_num_icode]]
                                    o.write(records_to_pdb_string(r) + '\n')
                                    break
                        else:
                            if r['resName'] == AA_Code[CDRs[numbering][res_num_icode]] or r['name'] in [' N  ', ' CA ',
                                                                                                        ' C  ', ' O  ']:
                                r['resName'] = AA_Code[CDRs[numbering][res_num_icode]]
                                o.write(records_to_pdb_string(r) + '\n')
    if not camelid:
        region_info = [('L', 'numbering_L', L_regions), ('H', 'numbering_H', H_regions)]
    else:
        region_info = [('H', 'numbering_H', H_regions)]


    for chain, numbering, regions in region_info:
        for R in regions:
            with open(prefix + R + '.pdb', 'w') as o:
                for line in open(prefix + '/template.' + R + '.pdb'):
                    r = map_pdb_string_to_records(line)
                    if r['type'] == "ATOM  ":
                        res_num = int(r['resSeq'])
                        res_num_icode = ('%s%s' % (res_num, r['iCode'])).strip()
                        if res_num_icode in CDRs[numbering] and chain == r['chainID']:
                            if regions[R][0] - 4 <= res_num <= regions[R][1] + 4:
                                if r['name'] in [' N  ', ' CA ', ' C  ', ' O  '] or r['resName'] == AA_Code[
                                    CDRs[numbering][res_num_icode]]:
                                    r['resName'] = AA_Code[CDRs[numbering][res_num_icode]]
                                    o.write(records_to_pdb_string(r) + '\n')

profit_templates = {'L': '''reference "%(prefix)s/template.light_heavy.pdb"
mobile "%(prefix)s/template.threaded.FRL.pdb"
atoms ca
zone L10-L23
zone L35-L49
zone L57-L66
zone L69-L88
zone L98-L100
fit
write "%(prefix)s/fitted.L.pdb"
quit''',
                    'H': '''reference "%(prefix)s/template.light_heavy.pdb"
mobile "%(prefix)s/template.threaded.FRH.pdb"
atoms ca
zone H10-H25
zone H36-H49
zone H66-H72
zone H74-H82
zone H83-H94
zone H103-H105
fit
write "%(prefix)s/fitted.H.pdb"
quit'''}

camelid_profit_template = {'H': '''reference "%(prefix)s/template.heavy.pdb"
mobile "%(prefix)s/template.threaded.FRH.pdb"
atoms ca
zone H10-H25
zone H36-H49
zone H66-H72
zone H74-H82
zone H83-H94
zone H103-H105
fit
write "%(prefix)s/fitted.H.pdb"
quit'''}


def superimpose_templates(CDRs, prefix):
    camelid = False
    if Options.camelid:
        camelid = True
    templates = ''
    if not camelid:
        templates = profit_templates
    else:
        templates = camelid_profit_template
    for chain in templates:
        f_name = prefix + 'profit-%s' % chain
        with open(f_name + '.in', 'w') as f:
            f.write(templates[chain] % vars())

        f_name = '\ '.join(f_name.split())
        commandline = '%s < %s.in > %s.out' % (Options.superimpose_profit, f_name, f_name)
        p = subprocess.Popen(commandline, shell=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
        output = p.communicate()[0].decode()
        res = p.poll()

        if res:
            print(commandline, output)
            sys.exit(1)

    pathPrefix = '\ '.join(vars()['prefix'].split())
    pathPrefix = pathPrefix[:-1] if pathPrefix.endswith('/') else pathPrefix
    commandline = ''
    if not camelid:
        commandline = 'cat {0}/fitted.L.pdb {0}/fitted.H.pdb > {0}/FR.pdb'.format(pathPrefix)
    else:
        commandline = 'cp {0}/fitted.H.pdb {0}/FR.pdb'.format(pathPrefix)
    p = subprocess.Popen(commandline, shell=True,
                         stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
    output = p.communicate()[0].decode()
    res = p.poll()

    if res:
        print(output)
        sys.exit(1)


def run_rosetta(CDRs, prefix, rosetta_bin, rosetta_platform, rosetta_database):
    addition_options = ''
    if Options.camelid:
        addition_options = ' -camelid'
    antibody_graft = rosetta_bin + '/antibody_graft.' + rosetta_platform
    if os.path.isfile(antibody_graft):
        print('\nRunning antibody_graft')
        # Sergey: Removing ' -restore_pre_talaris_2013_behavior' + \  because it lead to segfault on mpi-intel build
        commandline = 'cd "%s/details" && "%s" -database %s -overwrite -s FR.pdb' % \
                      (os.path.dirname(prefix), antibody_graft, rosetta_database) + \
                      ' -antibody::h3_no_stem_graft' + addition_options + ' -scorefile score-graft.sf'

        if Options.constant_seed:
            commandline += ' -run:constant_seed'
        if Options.quick:
            commandline += ' -run:benchmark -antibody:stem_optimize false'
        p = subprocess.Popen(commandline, shell=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
        output = p.communicate()[0].decode()
        res = p.poll()

        if Options.verbose or res:
            print(commandline, output)
        if res:
            print('Rosetta run terminated with Error!')
            sys.exit(1)
        model_file_prefix = 'grafted'
        shutil.move(prefix + 'details/FR_0001.pdb', prefix + model_file_prefix + '.pdb')
    else:
        print('Rosetta executable %s was not found, skipping Rosetta run...' % antibody_graft)
        return

    if Options.idealize:
        idealize_jd2 = rosetta_bin + '/idealize_jd2.' + rosetta_platform
        if os.path.isfile(idealize_jd2):
            print('Running idealize_jd2')
            commandline = 'cd "%s" && "%s" -database %s -overwrite' % (
            os.path.dirname(prefix), idealize_jd2, rosetta_database) + \
                          ' -fast -s %s.pdb -ignore_unrecognized_res -scorefile score-idealize.sf' % model_file_prefix
            p = subprocess.Popen(commandline, shell=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
            output = p.communicate()[0].decode()
            res = p.poll()

            if Options.verbose or res:
                print(commandline, output)
            if res:
                print('Rosetta run terminated with Error!  Commandline: %s' % commandline)
                sys.exit(1)
            shutil.move(prefix + model_file_prefix + '_0001.pdb', prefix + model_file_prefix + '.idealized.pdb')
            model_file_prefix += '.idealized'
        else:
            print('Rosetta executable %s was not found, skipping Rosetta run...' % idealize_jd2)
            return

    if Options.relax:
        relax = rosetta_bin + '/relax.' + rosetta_platform
        if os.path.isfile(relax):
            print('Running relax with all-atom constraint')
            commandline = 'cd "%s" && %s "%s" -database %s -overwrite' % (
            os.path.dirname(prefix), 'ulimit -t %s &&' % Options.timeout if Options.timeout else '', relax, rosetta_database) + \
                          ' -s %s.pdb -ignore_unrecognized_res -relax:fast -relax:constrain_relax_to_start_coords' % model_file_prefix + \
                          ' -relax:coord_constrain_sidechains -relax:ramp_constraints false -ex1 -ex2 -use_input_sc -scorefile score-relax.sf'
            p = subprocess.Popen(commandline, shell=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
            output = p.communicate()[0].decode()
            res = p.poll()

            if Options.verbose or res:
                print(commandline, output)
            if res:
                print('Rosetta run terminated with Error!  Commandline: %s' % commandline)
                sys.exit(1)
            shutil.move(prefix + model_file_prefix + '_0001.pdb', prefix + model_file_prefix + '.relaxed.pdb')
            model_file_prefix += '.relaxed'
        else:
            print('Rosetta executable %s was not found, skipping Rosetta run...' % relax)
            return

    shutil.copy(prefix + model_file_prefix + '.pdb', prefix + 'model.pdb')
    cter = kink_or_extend(CDRs)
    output_cter_constraint(cter, prefix,CDRs)


def kink_or_extend(CDRs):
    camelid = False
    if Options.camelid:
        camelid = True
    """Daisuke's H3 kink rules [citation]"""
    print('\nPreparing cter_constraint file for H3 modeling')
    if not camelid:
        L36 = AA_Code[CDRs['numbering_L']['36']]
        L46 = AA_Code[CDRs['numbering_L']['46']]
        L49 = AA_Code[CDRs['numbering_L']['49']]
    H93 = AA_Code[CDRs['numbering_H']['93']]
    H94 = AA_Code[CDRs['numbering_H']['94']]
    H99 = AA_Code[CDRs['H3'][-4:-3]]  # n-3
    H100 = AA_Code[CDRs['H3'][-3:-2]]  # n-2
    H101 = AA_Code[CDRs['H3'][-2:-1]]  # n-1
    len_h3 = len(CDRs['H3'])
    if not camelid:
        print('H3:', len_h3, CDRs['H3'], '\nKey residues: ', L36, L46, L49, H93, H94, H99, H100, H101)
    else:
        print('H3:', len_h3, CDRs['H3'], '\nKey residues: ', H93, H94, H99, H100, H101)

    # H3-rules
    if (H93 == 'ARG' or H93 == 'LYS') and (H94 == 'ARG' or H94 == 'LYS') and H101 == 'ASP':
        base = 'KINK'
    elif (H94 == 'ARG' or H94 == 'LYS') and H101 == 'ASP':
        if not camelid and (L46 == 'ARG' or L46 == 'LYS') and L36 != 'TYR':
            base = 'EXTEND'
        else:
            base = 'KINK'
    elif (H93 == 'ARG' or H93 == 'LYS') and H101 == 'ASP':
        base = 'KINK'
    elif H101 == 'ASP':
        if not camelid and (L49 == 'ARG' or L49 == 'LYS'):
            base = 'KINK'
        elif (H100 == 'MET' or H100 == 'PHE') and (H99 == 'ALA' or H99 == 'GLY'):
            base = 'KINK'
        else:
            base = 'EXTEND'
    elif (H100 == 'MET' or H100 == 'PHE') and (H99 == 'ALA' or H99 == 'GLY'):
        base = 'KINK'
    elif H100 == 'ASP' or H100 == 'ASN' or H100 == 'LYS' or H100 == 'ARG':
        base = 'EXTEND'
    elif len_h3 == 7:
        base = 'EXTEND'
    else:
        base = 'KINK'
    print('Predicted base conformation is', base)
    return base


# Dihedral CA 220 CA 221 CA 222 CA 223 SQUARE_WELL2 0.523 0.698 200; KINK
# Dihedral CA 220 CA 221 CA 222 CA 223 SQUARE_WELL2 2.704 0.523 100; EXTEND
def output_cter_constraint(base, prefix, CDRs):
    # Jianqing's original
    cnt = 0
    f = open(prefix + 'cter_constraint', 'w')
    for line in open(prefix + '/model.pdb'):
        if line[0:4] == 'ATOM':
            chain = line[21:22]
            atom = line[13:15]
            res_num = line[22:27].strip()

            if atom == 'CA':
                cnt += 1
            if chain == 'H' and atom == 'CA' and res_num == str(103):
                n1 = cnt - 3  # n-2 (H100X)
                n2 = cnt - 2  # n-1 (H101)
                n3 = cnt - 1  # n   (H102)
                n4 = cnt  # n+1 (H103)

                if base == 'KINK':
                    f.write('Dihedral CA ' + str(n1) + ' CA ' + str(n2) + ' CA ' + str(n3) + ' CA ' + str(
                        n4) + ' SQUARE_WELL2 0.523 0.698 200')
                elif base == 'EXTEND':
                    f.write('Dihedral CA ' + str(n1) + ' CA ' + str(n2) + ' CA ' + str(n3) + ' CA ' + str(
                        n4) + ' SQUARE_WELL2 2.704 0.523 100')
    f.close()

    # new python script
    if Options.kink_constraints:
        try: res = VhCdr3KinkConstraints.generate_constraints('%s/model.pdb' % prefix, CDRs['H3'])
        except Exception as e:
            print('ERROR making constraint file!', e)
            sys.exit(1)
#        commandline = '%s/kink_constraints.py %s/model.pdb' % (script_dir, prefix)
#        p = subprocess.Popen(commandline, shell=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
#        output = p.communicate()[0].decode()
#        res = p.poll()
#
#        if Options.verbose and not res:
#            print(commandline + '\n', output)
#        if res:
#            print(commandline + '\n', 'ERROR making constraint file: %s\n%s' % (res, output))
#            sys.exit(1)

class VhCdr3KinkConstraints(object):
    @staticmethod
    def find_vh_cdr3_location(pdbfilename,vh_cdr3):
        commandline = 'cat %s | awk \'/ATOM/ && $3 == "CA" {print $4}\' | tr \'\n\' \' \' | ' % pdbfilename  + \
        'sed \'s/ALA/A/g;s/CYS/C/g;s/ASP/D/g;s/GLU/E/g;s/PHE/F/g;s/GLY/G/g;s/HIS/H/g;s/ILE/I/g;s/LYS/K/g;' + \
        's/LEU/L/g;s/MET/M/g;s/ASN/N/g;s/PRO/P/g;s/GLN/Q/g;s/ARG/R/g;s/SER/S/g;s/THR/T/g;s/VAL/V/g;s/TRP/W/g;s/TYR/Y/g\'' + \
        ' | sed \'s/ //g\''
        p = subprocess.Popen(commandline, shell=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
        output = p.communicate()[0].decode()
        res = p.poll()
        if res: raise Exception('Wrong pdb file!')
        cdr3_begin = output.find(vh_cdr3) + 1
        if cdr3_begin == -1: raise Exception('Can\'t locate VH CDR3')
        cdr3_end = cdr3_begin + len(vh_cdr3) - 1
        return (cdr3_begin, cdr3_end)

    @staticmethod
    def kink_constraints(cdr3_begin, cdr3_end, outf):
        kink_begin = cdr3_end - 2
        kink_anion_res = cdr3_end - 1
        kink_kation_res = cdr3_begin - 1
        # CTER kink dihedral angle q
        outf.write("Dihedral CA %i CA %i CA %i CA %i SQUARE_WELL2 0.523 0.698 600\n" %
                (kink_begin, kink_begin + 1, kink_begin + 2, kink_begin + 3) )
        # CTER kink q bond distance
        #outf.write("AtomPair CA %i CA %i FLAT_HARMONIC 7.125 0.5 0.625\n" % (kr0,kr3) )
        # KD Hbond
        outf.write("AtomPair N %i O %i FLAT_HARMONIC 2.0 2.0 2.0\n" %
                (kink_anion_res, kink_kation_res) )
    

    @staticmethod
    def generate_constraints(pdbfilename, vh_cdr3):
        (cdr3_begin, cdr3_end) = VhCdr3KinkConstraints.find_vh_cdr3_location(pdbfilename, vh_cdr3)
        outfname = pdbfilename[:-4]+'.constr'
        outf = open(outfname, 'w')
        constraint = VhCdr3KinkConstraints.kink_constraints(cdr3_begin,cdr3_end,outf)


# Various filter function
def filter_by_sequence_homolog(k, results, cdr_query, cdr_info):
    camelid = False
    if Options.camelid:
        camelid = True
    if Options.verbose:
        print('filtering by sequence identity...')
    # print(results)
    if camelid:
        cdrs_names = ['H1', 'H2', 'H3']
    else:
        cdrs_names = ['L1', 'L2', 'L3', 'H1', 'H2', 'H3']

    for r in results[:]:
        pdb = r['subject-id']
        if k in cdrs_names and sid_checker(cdr_query[k], cdr_info[pdb][k]) >= sid_cutoff_cdr:
            results.remove(r)
            if Options.verbose and r not in results:
                print('Filter sequence_identity (%s%% cut-off), removing:%s %s_query:%s %s_info:%s %s%% identity' % \
                      (sid_cutoff_cdr, pdb, k, len(cdr_query[k]), k, len(cdr_info[pdb][k]),
            round(sid_checker(cdr_query[k], cdr_info[pdb][k]), 2)))

        elif not camelid and k in ['FRL'] and sid_checker(cdr_query[k], frl_info[pdb][k]) >= sid_cutoff_fr:
            results.remove(r)
            if Options.verbose and r not in results:
                print('Filter sequence_identity (%s%% cut-off), removing:%s %s_query:%s %s_info:%s %s%% identity' % \
                      (sid_cutoff_fr, pdb, k, len(cdr_query[k]), k, len(frl_info[pdb][k]),
            round(sid_checker(cdr_query[k], frl_info[pdb][k]), 2)))

        elif k in ['FRH'] and sid_checker(cdr_query[k], frh_info[pdb][k]) >= sid_cutoff_fr:
            results.remove(r)
            if Options.verbose and r not in results:
                print('Filter sequence_identity (%s%% cut-off), removing:%s %s_query:%s %s_info:%s %s%% identity' % \
                      (sid_cutoff_fr, pdb, k, len(cdr_query[k]), k, len(frh_info[pdb][k]),
            round(sid_checker(cdr_query[k], frh_info[pdb][k]), 2)))

        elif not camelid and k in ['light_heavy'] and float(r['%-identity']) >= sid_cutoff_fr - 5:
            results.remove(r)
            if Options.verbose and r not in results:
                print('Filter sequence_identity (%s%% cut-off), removing:%s %s %s%% identity' % \
                      (sid_cutoff_fr, pdb, k, r['%-identity']))
        elif camelid and k in ['heavy'] and float(r['%-identity']) >= sid_cutoff_fr - 5:
            results.remove(r)
            if Options.verbose and r not in results:
                print('Filter sequence_identity (%s%% cut-off), removing:%s %s %s%% identity' % \
                      (sid_cutoff_fr, pdb, k, r['%-identity']))


def sid_checker(seq_q, seq_t):
    seq_lenq = len(seq_q)
    seq_lent = len(seq_t)

    if seq_lenq == seq_lent:
        num_res = 0
        for i in range(0, seq_lenq):
            if seq_q[i] == seq_t[i]:
                num_res += 1
        ratio = float(num_res) / seq_lenq * 100
    else:
        ratio = 0

    return ratio


def filter_by_sequence_length(k, results, cdr_query, cdr_info):
    if Options.verbose:
        print('filtering by sequence length...')

    for r in results[:]:
        pdb = r['subject-id']
        if k in ['heavy', 'light_heavy']:
            return
        elif k in ['L1', 'L2', 'L3', 'H1', 'H2', 'H3'] and len(cdr_query[k]) != len(cdr_info[pdb][k]):
            results.remove(r)
            if Options.verbose and r not in results:
                print('Filter sequence_length, removing:%s %s_query:%s %s_info:%s' % (
            pdb, k, len(cdr_query[k]), k, len(cdr_info[pdb][k])))
        elif k == 'light' and 'light_lenght' in cdr_info[pdb] and len(cdr_query[k]) != int(
                cdr_info[pdb]['light_lenght']):
            results.remove(r)
        elif k == 'FRH':
            template_length = 61 if pdb == 'pdb2x7l_chothia.pdb' else 63
            if not len(cdr_query[k]) == template_length:
                results.remove(r)
            if Options.verbose and r not in results:
                print('Filter sequence_length, removing:%s %s_query:%s %s_info:%s' % (
            pdb, k, len(cdr_query[k]), k, template_length))
            # elif k == 'FRH'  and  (not  len(cdr_query[k])-8 <= 67 <= len(cdr_query[k])+8): results.remove(r)
        elif k == 'FRL':
            template_length = 60 if pdb == 'pdb3h0t_chothia.pdb' else 58
            if not len(cdr_query[k]) == template_length:
                results.remove(r)
            if Options.verbose and r not in results:
                print('Filter sequence_length, removing:%s %s_query:%s %s_info:%s' % (
            pdb, k, len(cdr_query[k]), k, template_length))
            # if not  len(cdr_query[k])-8 <= template_length <= len(cdr_query[k])+8: results.remove(r)


def filter_by_alignment_length(k, results, cdr_query, cdr_info):
    for r in results[:]:
        pdb = r['subject-id']
        if k == 'H3' and int(r['alignment-length']) < 0.10 * len(cdr_info[pdb]['H3']):
            results.remove(r)
        if k == 'H2' and int(r['alignment-length']) < 0.55 * len(cdr_info[pdb]['H2']):
            results.remove(r)
        if k in ['L1', 'L2', 'L3', 'H1'] and int(r['alignment-length']) < 0.70 * len(cdr_info[pdb][k]):
            results.remove(r)
        if Options.verbose and r not in results:
            print('Filter alignment_length removing:%s %s_query:%s alignment-length:%s ' % \
                  (pdb, k, len(cdr_info[pdb][k]), r['alignment-length']))


def filter_by_template_resolution(k, results, cdr_query, cdr_info):
    for r in results[:]:
        pdb = r['subject-id']
        if float(cdr_info[pdb]['resolution']) > 2.8:
            results.remove(r)
        if Options.verbose and r not in results:
            print('Filter template_resolution, removing:%s resolution:%s' % (
        pdb, cdr_info[pdb]['resolution']))


def filter_by_template_bfactor(k, results, cdr_query, cdr_info):
    bfactor = {}

    if k in ['L1', 'L2', 'L3', 'H1', 'H2', 'H3']:
        for line in open(_script_path_ + '/info/list_bfactor50'):
            for i, e in enumerate(['L1', 'L2', 'L3', 'H1', 'H2', 'H3']):
                if k == e:
                    bfactor['pdb' + line.split()[0] + '_chothia.pdb'] = line.split()[i + 1] == 'True'

        for r in results[:]:
            pdb = r['subject-id']
            if bfactor.get(pdb, False):
                results.remove(r)
            if Options.verbose and r not in results:
                print('Filter B-factor50, removing:%s' % pdb)


def filter_by_outlier(k, results, cdr_query, cdr_info):
    outlier = {}
    for line in open(_script_path_ + '/info/outlier_list'):
        outlier[tuple(line.split()[:2])] = line.split()[2] == 'true'
    for r in results[:]:
        pdb = r['subject-id']

        if outlier.get((pdb, k), False):
            results.remove(r)
        if Options.verbose and r not in results:
            print('Filter outlier, removing:%s' % pdb)


AA_Code = dict(A='ALA', V='VAL', L='LEU', I='ILE', P='PRO', W='TRP', F='PHE', M='MET', G='GLY', S='SER', T='THR',
               Y='TYR', C='CYS', N='ASN', Q='GLN', K='LYS', R='ARG', H='HIS', D='ASP', E='GLU')

# PDB Reader code
# PDB record maps, only ATOMS here for now
PDB_Records = {
    "ATOM  ": {
        "type": (1, 6),
        "serial": (7, 11),  # Integer
        "name": (13, 16),  # Atom
        "altLoc": (17, 17),  # Character
        "resName": (18, 20),  # Residue name
        "chainID": (22, 22),  # Character
        "resSeq": (23, 26),  # Integer
        "iCode": (27, 27),  # AChar
        "x": (31, 38),  # Real(8.3)
        "y": (39, 46),  # Real(8.3)
        "z": (47, 54),  # Real(8.3)
        "occupancy": (55, 60),  # Real(6.2)
        "tempFactor": (61, 66),  # Real(6.2)  #///"segID",     Field(73, 76),
        "element": (77, 78),  # LString(2)
        "charge": (79, 80)  # LString(2)
    },
    "UNKNOWN": {
        "type": (1, 6),
        "info": (7, 80)
    }
}


def map_pdb_string_to_records(s):
    rtype = s[:6]
    F = PDB_Records.get(rtype, PDB_Records["UNKNOWN"])
    R = {}
    for f in F:
        R[f] = s[F[f][0] - 1: F[f][1]]
    return R


def records_to_pdb_string(records):
    res = [' '] * 80
    for f in records:
        res[PDB_Records['ATOM  '][f][0] - 1: PDB_Records['ATOM  '][f][1]] = records[f]
    return ''.join(res)


def self_test():
    if os.path.isdir(Options.self_test_dir):
        print('Removing old self-test-dir %s...' % Options.self_test_dir)
        shutil.rmtree(Options.self_test_dir)
    os.makedirs(Options.self_test_dir)  # if not os.path.isdir(Options.self_test_dir):

    tests = [d for d in os.listdir('test') if d != '.svn']
    print('Checking [%s] targets: %s...' % (len(tests), tests))
    for t in tests:
        test_dir = Options.self_test_dir + t + '/'
        os.makedirs(test_dir)
        if Options.verbose:
            print('Testing target: %s...' % t)
        light_chain = read_fasta_file('test/%s/query_l.fasta' % t)
        heavy_chain = read_fasta_file('test/%s/query_h.fasta' % t)
        answers = json.load(open('test/%s/%s.json' % (t, t)))
        answers['numbering_L'] = open('test/%s/numbering_L.txt' % t).read()
        answers['numbering_H'] = open('test/%s/numbering_H.txt' % t).read()

        if Options.verbose:
            print('light_chain:', light_chain)
            print('heavy_chain:', heavy_chain)

        CDRs = IdentifyCDRs(light_chain, heavy_chain)
        if Options.verbose: print('CDR:', json.dumps(CDRs, sort_keys=True, indent=2))
        CDRs.update(Extract_FR_CDR_Sequences(**CDRs))

        write_results(CDRs, prefix=test_dir)
        CDRs['numbering_L'] = open(test_dir + '/numbering_L.txt').read()
        CDRs['numbering_H'] = open(test_dir + '/numbering_H.txt').read()

        for a in sorted(answers.keys()):  # ['L1', 'L2', 'L3', 'H1', 'H2', 'H3']
            if answers[a] != CDRs[a]:
                print('ERROR: target=%s field %s is not equal!!!\nexpected:%s\n    got:%s' % (t, a, answers[a], CDRs[a]))
                sys.exit(1)
            elif Options.verbose:
                print('  %s: OK' % a)

        if Options.verbose:
            print('Testing target: %s... OK\n' % t)

    print('All tests passed!')
    return


if __name__ == "__main__":
    main(sys.argv)
