#!/usr/bin/env python
# encoding: utf-8
'''
bigDump -- Do a large Entry query using BioPython

bigDump performs a large fetch and export from NCBI using a list of gi numbers

@author: Aaron Kitzmiller <aaron_kitzmiller@harvard.edu>
@copyright: 2016 The Presidents and Fellows of Harvard College. All rights reserved.
@license: GPL v2.0
'''

import sys, os, time, traceback

from Bio import Entrez
from urllib2 import HTTPError
import logging


from argparse import ArgumentParser
from argparse import RawDescriptionHelpFormatter

__all__ = []
__version__ = 0.1
__date__ = '2016-11-03'
__updated__ = '2016-11-03'

DEBUG = 1
PROFILE = 0


def main(argv=None): # IGNORE:C0111
    '''Command line options.'''

    if argv is None:
        argv = sys.argv
    else:
        sys.argv.extend(argv)

    program_name = os.path.basename(sys.argv[0])
    program_version = "v%s" % __version__
    program_build_date = str(__updated__)
    program_version_message = '%%(prog)s %s (%s)' % (program_version, program_build_date)
    program_shortdesc = __import__('__main__').__doc__.split("\n")[1]
    program_license = '''%s

  Created by Aaron Kitzmiller on %s.
  Copyright 2016 The Presidents and Fellows of Harvard College. All rights reserved.

  Licensed under the Gnu General Public License v2.0

  Distributed on an "AS IS" basis without warranties
  or conditions of any kind, either express or implied.

USAGE
''' % (program_shortdesc, str(__date__))

    try:
        # Setup argument parser
        parser = ArgumentParser(description=program_license, formatter_class=RawDescriptionHelpFormatter)
        parser.add_argument("-v", "--verbose", dest="verbose", action="store_true", help="set verbosity level [default: %(default)s]")
        parser.add_argument('-V', '--version', action='version', version=program_version_message)
        parser.add_argument("gis",metavar="GIFILE",help="File containing a newline-separated list of gi numbers")
        parser.add_argument('-e','--email',required=True,help="Email address for NCBI Entrez")
        parser.add_argument('-d','--database',help="NCBI database. Default 'protein'.")
        parser.set_defaults(database='protein',verbose=False)

        # Process arguments
        args = parser.parse_args()

        verbose = args.verbose 
        if verbose:
            logging.basicConfig(format='%(message)s',level=logging.INFO)
        else:
            logging.basicConfig(format='%(message)s',level=logging.ERROR)

                        
        if not os.path.exists(args.gis):
            raise Exception("gi list file %s does not exist" % args.gis)
        
        db = args.database
        
        Entrez.email = args.email
        
        # Read the gis
        logging.info('Reading gis')       
        gis = []
        with open(args.gis,'r') as f:
            for line in f:
                gis.append(line.strip())
        logging.info('Read %d gis from %s' % (len(gis),args.gis))
        # Post to Entrez
        logging.info('POSTing to Entrez')            
        result = Entrez.read(Entrez.epost('protein',id=','.join(gis)))
        querykey = result['QueryKey']
        webenv = result['WebEnv']
        logging.info('Done POSTing to Entrez')
            
        # Download a batch at a time
        count = len(gis)
        batchsize = 200
        for start in range(0, count, batchsize):
            end = min(count, start + batchsize)
            logging.info("Going to download record %i to %i" % (start+1, end))
            attempt = 1
            while attempt <= 3:
                try:
                    logging.info('Attempt %d' % attempt)
                    handle = Entrez.efetch(db=db, retmode="xml",
                        retstart=start, retmax=batchsize,
                        webenv=webenv, query_key=querykey)
                    attempt = 4
                except HTTPError as err:
                    if 500 <= err.code <= 599:
                        logging.info("Received error from server %s" % err)
                        logging.info("Attempt %i of 3" % attempt)
                        attempt += 1
                        time.sleep(15)
                    else:
                        raise
            for r in Entrez.parse(handle):
                try:
                    gi=int([x for x in r['GBSeq_other-seqids'] if "gi" in x][0].split("|")[1])
                except ValueError:
                    gi=None
                print '>gi|{gi}|gb|{acc}| {definition}\n{seq}'.format(gi=gi,acc=r["GBSeq_primary-accession"],definition=r["GBSeq_definition"],seq=r["GBSeq_sequence"].upper())

            handle.close()

        return 0

    except KeyboardInterrupt:
        ### handle keyboard interrupt ###
        return 0
    
    except Exception as e:
        sys.stderr.write('Error %s\n%s\n' % (str(e),traceback.format_exc()))
        return 2

if __name__ == "__main__":
    if DEBUG:
        sys.argv.append("-v")
    if PROFILE:
        import cProfile
        import pstats
        profile_filename = 'bigDump_profile.txt'
        cProfile.run('main()', profile_filename)
        statsfile = open("profile_stats.txt", "wb")
        p = pstats.Stats(profile_filename, stream=statsfile)
        stats = p.strip_dirs().sort_stats('cumulative')
        stats.print_stats()
        statsfile.close()
        sys.exit(0)
    sys.exit(main())
