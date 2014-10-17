#!/usr/bin/env python
import os
import sys
from optparse import OptionParser

if __name__=="__main__":
    """Create a new directory in the filetree.

    This script creates the directory outputdir.  It takes care of adding
    intermediate directories as needed.  This script doesn't delete anything
    already present in the filesytem.
    """

    parser = OptionParser()
    parser.add_option("-o", "--outputdir", dest="outputdir",
                      help="output directory to be created",
                      default='output', metavar="OUTPUTDIR")
    # TODO - maybe add in somethign that prints a friendly help message or the
    # docstring?

    (options, args) = parser.parse_args()

    # create the outputdirectory
    dir = options.outputdir
    if not os.path.exists( dir ):
        os.makedirs( dir )

    # TODO - copy over parameters and other files used to run this example
