#!/usr/bin/env python
# encoding: utf-8
'''
bigDump -- Do a large Entry query using BioPython

bigDump performs a large fetch and export from NCBI using a list of gi numbers

@author: Aaron Kitzmiller <aaron_kitzmiller@harvard.edu>
@copyright: 2016 The Presidents and Fellows of Harvard College. All rights reserved.
@license: GPL v2.0
'''

import sys, os, time

from Bio import Entrez
from urllib2 import HTTPError
import logging

logger = logging.getLogger()
logger.setLevel(logging.ERROR)


from argparse import ArgumentParser
from argparse import RawDescriptionHelpFormatter

__all__ = []
__version__ = 0.1
__date__ = '2016-11-03'
__updated__ = '2016-11-03'

DEBUG = 1
PROFILE = 0


BATCH_SIZE = 200 # Size of the download batch


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
        parser.set_defaults(database='protein')

        # Process arguments
        args = parser.parse_args()

        verbose = args.verbose 
        if verbose:
            logger.setLevel(logging.INFO)
                        
        if not os.path.exists(args.gis):
            raise Exception("gi list file %s does not exist" % args.gis)
        
        db = args.database
        
        Entrez.email = args.email
        
        # Read the gis
        logger.info('Reading gis')       
        gis = []
        with open(args.gis,'r') as f:
            gis = f.read()    
        logger.info('Read %d gis from %s' % (len(gis),args.gis))
        
        # Post to Entrez
        logger.info('POSTing to Entrez')            
        result = Entrez.ePost('protein',id=','.join(gis))
        querykey = result['QueryKey']
        webenv = result['WebEnv']
        logger.info('Done POSTing to Entrez')
            
        # Download a batch at a time
        count = len(gis)
        for start in range(0, count, BATCH_SIZE):
            end = min(count, start + BATCH_SIZE)
            logger.info("Going to download record %i to %i" % (start+1, end))
            attempt = 1
            while attempt <= 3:
                try:
                    handle = Entrez.efetch(db=db, rettype="fasta", retmode="text",
                        retstart=start, retmax=BATCH_SIZE,
                        webenv=webenv, query_key=querykey)
                except HTTPError as err:
                    if 500 <= err.code <= 599:
                        logger.info("Received error from server %s" % err)
                        logger.info("Attempt %i of 3" % attempt)
                        attempt += 1
                        time.sleep(15)
                    else:
                        raise
            data = handle.read()
            handle.close()
            print data

        return 0

    except KeyboardInterrupt:
        ### handle keyboard interrupt ###
        return 0
    
    except Exception, e:
        if DEBUG:
            raise(e)
        indent = len(program_name) * " "
        sys.stderr.write(program_name + ": " + repr(e) + "\n")
        sys.stderr.write(indent + "  for help use --help")
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