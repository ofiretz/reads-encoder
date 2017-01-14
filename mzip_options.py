from optparse import OptionParser
from mzpy_options import *
from mzpy_math import ispower2
import vlc
import re

def parse_code_type(option, opt_str, value, parser):
    assert value is None
    value = ["", []]

    try:
        code = parser.rargs[0]
    except IndexError:
        raise CLIError("missing argument to '{0}'".format(opt_str))


    if code not in compression_methods:
        raise CLIError("unknown argument {0} to option {1}".format(code
            , opt_str))

    value[0] = code

    if code == "golomb":
        try:
            m = parser.rargs[1]
        except IndexError:
            raise CLIError("missing argument to '{0} {1}'".format(opt_str
                , code))
        if not re.match(r'^\d+$', m):
            raise CLIError("argument to '{0} {1}' must be integer"\
                .format(opt_str, code))
        m = int(m)
        value[1].append(m)

    num_read = 1 + len(value[1])

    del parser.rargs[:num_read]
    setattr(parser.values, option.dest, value)

def parse_commandline(argv):
    cl_parser = OptionParser()

    cl_parser.add_option("-f", "--file", type="string"
        , dest="infile", help="input file")

    cl_parser.add_option("-r", "--rpos", dest="rpos"
        , action="callback", callback=parse_code_type
        , help="encoding of read positions") 
    
    cl_parser.add_option("-v", "--vpos", dest="vpos"
        , action="callback", callback=parse_code_type
        , help="encoding of variation positions")

    cl_parser.add_option("-p", "--paired", action="store_true", dest="paired"
        , default=False, help="indicate paired end reads")

    cl_parser.add_option("-a", "--paircode", dest="poff"
        , default=['gamma', []]
        , action="callback", callback=parse_code_type
        , help="encoding of read mates offset")

    if len(argv) == 1:
        cl_parser.print_help()
        raise CLIEmpty()

    options, args = cl_parser.parse_args()

    if not options.infile:
        raise CLIError("option '--file' is mandatory")

    return options, args
