#! /usr/bin/env python

# {{{ imports
from __future__ import print_function, division
from sys import argv, stderr
from bitarray import bitarray
from functools import partial
import numpy as np
from collections import defaultdict
from operator import itemgetter
import Huffman3 as huff
from math import ceil
import csv
import os

from mzpy_globals import *
from vlc import *
from mzpy_options import *
from mzip_options import *
from bases import *
# }}}

try:
    options, args = parse_commandline(argv)
except CLIEmpty:
    exit(2)
except CLIError as e:
    print("{0}: {1}".format(argv[0], e.value))
    exit(2)

# {{{ BookKeeper
class BookKeeper(object):
    def __init__(self):
        self.misc = 0
        self.rlens = 0
        self.rpos = 0
        self.strands = 0
        self.matchflags = 0
        self.spos = 0
        self.sbases = 0
        self.sflags = 0
        self.ipos = 0
        self.ibases = 0
        self.iflags = 0
        self.dpos = 0
        self.dlens = 0
        self.dflags = 0
        self.pflags = 0
        self.poffs = 0
        self.quals = 0
    def show(self):
        print("***** storage in bits: *****")
        for type, nbits in self.__dict__.items():
            print(type, nbits)
        tot = sum(self.__dict__.values())
        tot_bytes = int(ceil(tot/8))
        print("TOTAL: {0} ({1} Bytes)".format(tot, tot_bytes))
        if tot_bytes != os.path.getsize(out_fn):
            print("error in book keeping, num bytes don't match", file=stderr)
# }}}
bookkeeping = BookKeeper()

class CommaCoder(object):
    def __init__(self, encoderfunc=None, commacode=None):
        self.encode = encoderfunc
        self.comma = commacode

# {{{ populate_bitwriter()
"""
rpos can be 0, but we don't need a comma code, i.e. gamma gets offset by 1
    and rice by 0

poff is always > 0, so no offset necessary

for the others we need a comma code and in some cases (first var position,
    joined indels) values can be 0, so we offset gamma by 2 and rice by 1,
    where gamma=1 / rice=0 will serve as comma
"""
def populate_bitwriter(opts):
    writer = defaultdict(CommaCoder)
    dispatch = {'gamma': gamma_encode,
                'golomb': golomb_encode, 
                'rice': rice_encode}
    for opt in compression_opts:
        method, params = getattr(opts, opt)
        if method == "gamma":
            off = 2 # default
            if opt == "rpos":
                off = 1
            elif opt == "poff":
                off = 0
            else:
                writer[opt].comma = dispatch['gamma'](1)
            writer[opt].encode = partial(dispatch['gamma'], off=off)
        if method == "golomb":
            m = params[0]
            if ispower2(m):
                func = dispatch["rice"]
            else:
                func = dispatch["golomb"]
            off = 1
            if opt in ["rpos", "poff"]:
                off = 0
            else:
                writer[opt].comma = func(0, m)
            writer[opt].encode = partial(func, m=m, off=off)
    return writer
# }}}

# {{{ myopen
def myopen(fname, mode="r"):
    try:
        fo = open(fname, mode)
    except IOError as err:
        print("{0}: open() failed -- {1}".format(argv[0], err), file=stderr)
        exit(2)
    return(fo)
# }}}

# {{{ write_magic()
def write_magic(stream):
    fourbytes = np.array([0], dtype=np.uint32)
    fourbytes[0] = MZ_MAGIC_32
    stream.write(fourbytes[0])
    bookkeeping.misc += 4 * 8
# }}}

# {{{ write_padlen()
def write_padlen(padding, stream):
    onebyte = np.array([0], dtype=np.uint8)
    onebyte[0] = len(padding)
    stream.write(onebyte[0])
    bookkeeping.misc += 8
# }}}

# {{{ write_main_flags()
def write_main_flags(stream
                     , opts
                     , same_length
                     , jumps_present
                     , identities_present
                     , qualities_present):
    twobytes = np.array([0], dtype=np.uint16)
    fourbytes = np.array([0], dtype=np.uint32)

    shift = 32

    val = 1 if same_length else 0 
    shift -= 1
    fourbytes[0] |= val << shift

    shift -= 1
    val = 1 if opts.paired else 0
    fourbytes[0] |= val << shift

    shift -= 1
    val = 1 if jumps_present else 0
    fourbytes[0] |= val << shift

    shift -= 1
    val = 1 if identities_present else 0
    fourbytes[0] |= val << shift

    shift -= 1
    val = 1 if qualities_present else 0
    fourbytes[0] |= val << shift

    numbits = 1 + int(floor(log(len(compression_methods)-1, 2)))
    for opt in compression_opts:
        method, params = getattr(opts, opt)
        val = compression_methods.index(method)
        shift -= numbits
        fourbytes[0] |= val << shift

    stream.write(fourbytes[0])
    bookkeeping.misc += 4 * 8

    # compression params:
    for opt in compression_opts:
        method, params = getattr(opts, opt)
        if method in ["golomb", "rice"]:
            m = params[0]
            twobytes[0] = m
            stream.write(twobytes[0])
            bookkeeping.misc += 2 * 8
# }}}

# {{{ write_seqlen()
def write_seqlen(stream, seqlen):
    fourbytes = np.array([0], dtype=np.uint32)
    fourbytes[0] = seqlen
    stream.write(fourbytes[0])
    bookkeeping.rlens += 4 * 8
# }}}

# {{{ write_len_num_pairs
def write_len_num_pairs(stream, len_num_pairs):
    twobytes = np.array([0], dtype=np.uint16)
    fourbytes = np.array([0], dtype=np.uint32)
    eightbytes = np.array([0], dtype=np.uint64)

    assert len(len_num_pairs) <= 2**32
    fourbytes[0] = len(len_num_pairs)
    stream.write(fourbytes[0])
    bookkeeping.rlens += 4 * 8

    for length, num in len_num_pairs:
        twobytes[0] = length
        eightbytes[0] = num
        stream.write(twobytes[0])
        stream.write(eightbytes[0])
        bookkeeping.rlens += (2 + 8) * 8
# }}}

# {{{ write_qual_num_pairs
def write_qual_num_pairs(stream, qual_num_pairs):
    onebyte = np.array([0], dtype=np.uint8)
    eightbytes = np.array([0], dtype=np.uint64)

    assert len(qual_num_pairs) <= 2**8
    onebyte[0] = len(qual_num_pairs)
    stream.write(onebyte[0])
    bookkeeping.quals += 8

    for qual, num in qual_num_pairs:
        onebyte[0] = ord(qual)
        eightbytes[0] = num
        stream.write(onebyte[0])
        stream.write(eightbytes[0])
        bookkeeping.quals += (1 + 8) * 8
# }}}

# {{{ write_jumps
def write_jumps(stream, jumps):
    fourbytes = np.array([0], dtype=np.uint32)
    eightbytes = np.array([0], dtype=np.uint64)

    assert len(jumps) <= 2**32
    fourbytes[0] = len(jumps)
    stream.write(fourbytes[0])
    bookkeeping.rpos += 4 * 8

    for rnum, rpos in jumps:
        eightbytes[0] = rnum
        fourbytes[0] = rpos
        stream.write(eightbytes[0])
        stream.write(fourbytes[0])
        bookkeeping.rpos += (4 + 8) * 8
# }}}

# {{{ get_median()
# input is dict with counts as values
def get_median(d):
    x_cnt = d.items()
    x_cnt.sort(key=itemgetter(0))
    tot_cnt = sum(map(itemgetter(1), x_cnt))
    cnt_accum = 0
    median = 0
    for x, cnt in x_cnt:
        cnt_accum += cnt
        median = x
        if cnt_accum > tot_cnt//2:
            break
    return median
# }}}

csvfile = open(options.infile, "rb")
reader = csv.reader(csvfile, delimiter="\t")

# {{{ get seq lengths, and position distributions, pair offsets
seq_lengths = defaultdict(int)
lastreadpos = 0
readposs = defaultdict(int)
varposs = defaultdict(int)
poffs = defaultdict(int)
pair_info = {}
identities_present = False
qual_cnts = defaultdict(int)
i = 0
for id, readpos, strand, readseq, varstr in reader:
    i += 1
    seqlen = len(readseq)
    seq_lengths[seqlen] += 1

    if not options.rpos:
        readpos = int(readpos)
        relreadpos = readpos - lastreadpos
        lastreadpos = readpos
        readposs[relreadpos] += 1

    if not options.vpos:
        lastvarpos = 0
        var = eval(varstr)
        for varpos, vartype, varmisc in var:
            varpos = int(varpos)
            relvarpos = varpos - lastvarpos
            lastvarpos = varpos
            varposs[relvarpos] += 1
            if vartype == "S" and varmisc[0] == varmisc[1]:
                identities_present = True
            #if quals != None:
             #   if vartype == "S":
              #      qual_cnts[quals] += 1
               # elif vartype == "I":
                #    for q in quals:
                 #       qual_cnts[q] += 1

    if options.paired:
        if id in pair_info:
            poff = i - pair_info[id]['i']
            pair_info[id] = {'ispaired': True, 'i': poff}
            poffs[poff] += 1
        else:
            pair_info[id] = {'ispaired': False, 'i': i}

if qual_cnts:
    qualities_present = True
else:
    qualities_present = False
# TODO find better estimator than median
if not options.rpos:
    median = get_median(readposs)
    if median < 2:
        median = 2
    options.rpos = ['golomb', [median]]

# TODO this was added last minute and is a bit clunky
# --> for big position differences, store (read_number, relative position)
# in 8+4 bytes instead of VLC
if options.rpos:
    method, params = options.rpos
    if method == "gamma":
        rpos_encoder = partial(gamma_encode, off=1)
    elif method == "golomb":
        m = params[0]
        if ispower2(m):
            rpos_encoder = partial(rice_encode, m=m)
        else:
            rpos_encoder = partial(golomb_encode, m=m)
else:
    m = median
    if ispower2(m):
        rpos_encoder = partial(rice_encode, m=m)
    else:
        rpos_encoder = partial(golomb_encode, m=m)
    
# need 12bytes=96bits to encode relative position (see above)
# so we'll use this for everything above 100 bits...
# FIXME shouldn't actually compute actual code here, but just length

csvfile.seek(0)
reader = csv.reader(csvfile, delimiter="\t")
i = 0
lastreadpos = 0
jumps = [] # pairs: (read_number, relpos)
for id, readpos, strand, readseq, varstr in reader:
    i += 1
    readpos = int(readpos)
    relreadpos = readpos - lastreadpos
    lastreadpos = readpos

    rpos_code = rpos_encoder(relreadpos)
    if len(rpos_code) >= 100:
        jumps.append((i, relreadpos))
if jumps:
    jumps_present = True
else:
    jumps_present = False

if not options.vpos:
    median = get_median(varposs)
    if median < 2:
        median = 2
    options.vpos = ['golomb', [median]]

if options.paired:
    median = get_median(poffs)
    if median < 2:
        median = 2
    options.poff = ['golomb', [median]]

del readposs
del varposs
del poffs
# }}}

len_num_pairs = seq_lengths.items()
#print(len_num_pairs)

if len(seq_lengths) == 1:
    same_length = True
    seqlen = len_num_pairs[0][0]
else:
    same_length = False
    tot = sum(map(itemgetter(1), len_num_pairs))
    probs = [(str(x[0]), x[1]/tot) for x in len_num_pairs]
    symbols = huff.makenodes(probs)
    root = huff.iterate(symbols)

if qualities_present:
    qual_num_pairs = qual_cnts.items()
    tot_q = sum( map(itemgetter(1), qual_num_pairs) )
    probs_q = [ (x[0], x[1]/tot_q) for x in qual_num_pairs]
    symbols_q = huff.makenodes(probs_q)
    root_q = huff.iterate(symbols_q)

bitwriter = populate_bitwriter(options)

out_fn = "{0}.mz".format(options.infile)
fout = myopen(out_fn, "wb")

#print("identities{0} present".format("" if identities_present else " not"))
#print("qualities{0} present".format("" if qualities_present else " not"))

write_magic(fout)
write_main_flags(fout
                 , options
                 , same_length
                 , jumps_present
                 , identities_present
                 , qualities_present)

buf = bitarray()

if same_length:
    write_seqlen(fout, seqlen)
else:
    write_len_num_pairs(fout, len_num_pairs)

if qualities_present:
    write_qual_num_pairs(fout, qual_num_pairs)

#print(len(jumps), " jumps present")

if jumps_present:
    write_jumps(fout, jumps)

csvfile.seek(0)
reader = csv.reader(csvfile, delimiter="\t")

lastreadpos = 0
i = 0
for id, readpos, strand, readseq, varstr in reader:
    i += 1
    readpos = int(readpos)
    relpos = readpos - lastreadpos
    assert relpos >= 0
    lastreadpos = readpos
    seqlen = len(readseq)
    var = eval(varstr)

    if not same_length:
        seqlen_ID = huff.encodesymbol(str(seqlen), symbols)
        buf.extend(seqlen_ID)
        bookkeeping.rlens += len(seqlen_ID)

    if jumps and jumps[0][0] == i:
        jumps.pop(0)
    else:
        poscode = bitwriter["rpos"].encode(relpos)
        buf.extend(poscode)
        bookkeeping.rpos += len(poscode)

    buf.extend(strand)
    bookkeeping.strands += 1

    if not var:
        buf.extend('1') # exact match
        bookkeeping.matchflags += 1
    else:
        buf.extend('0')
        bookkeeping.matchflags += 1
        lastvarpos = 0
        for varpos, vartype, varmisc in var:
            relvarpos = varpos - lastvarpos
            assert relvarpos >= 0
            lastvarpos = varpos

            poscode = bitwriter["vpos"].encode(relvarpos)            
            buf.extend(poscode)

            if vartype == "S":
                bookkeeping.spos += len(poscode)
                buf.extend('0')
                bookkeeping.sflags += 1
                read_base, ref_base = varmisc
                if identities_present:
                    if read_base == ref_base:
                        buf.extend("1")
                        bookkeeping.sbases += 1
                    else:
                        buf.extend("0")
                        subcode = ((base_codes[read_base] - base_codes[ref_base]) % 5) - 1
                        buf.extend("{0:02b}".format(subcode))
                        bookkeeping.sbases += 3
                else:
                    subcode = ( (base_codes[read_base] - base_codes[ref_base]) % 5 ) - 1
                    buf.extend("{0:02b}".format(subcode))
                    bookkeeping.sbases += 2
                if qualities_present:
                    qual_code = huff.encodesymbol(quals, symbols_q)
                    buf.extend(qual_code)
                    bookkeeping.quals += len(qual_code)

            if vartype == "I":
                bookkeeping.ipos += len(poscode)
                buf.extend('10')
                bookkeeping.iflags += 2
                bases_code = ''
                for base in varmisc:
                    bases_code += base2tbe[base]
                bases_code += base2tbe['eof']
                buf.extend(bases_code)
                bookkeeping.ibases += len(bases_code)
                if qualities_present:
                    for q in quals:
                        qual_code = huff.encodesymbol(q, symbols_q)
                        buf.extend(qual_code)
                        bookkeeping.quals += len(qual_code)

            if vartype == "D":
                bookkeeping.dpos += len(poscode)
                buf.extend('11')
                bookkeeping.dflags += 2
                len_code = gamma_encode(varmisc)
                buf.extend(len_code)
                bookkeeping.dlens += len(len_code)

        comma = bitwriter["vpos"].comma
        buf.extend(comma)
        bookkeeping.misc += len(comma)
               
    # FIXME: all flags are dummy values; BAM paired ends not supported yet     
    if options.paired:
        if pair_info[id]['ispaired']:
            if 'mate_processed' not in pair_info[id]:
                pInfo = '111' # is_paired, orientation, sequenced strand
                buf.extend(pInfo)
                poff = pair_info[id]['i']
                poff_code = bitwriter["poff"].encode(poff)
                buf.extend(poff_code)
                bookkeeping.pflags += 3
                bookkeeping.poffs += len(poff_code)                
                pair_info[id]['mate_processed'] = True
            else:
                pInfo = '1' # just sequenced strand
                buf.extend(pInfo)
                bookkeeping.pflags += 1

csvfile.close()
if buf.length() % 8 == 0:
    padding = ''
else:
    padding = '0' *  ( 8 - ( buf.length() % 8 ) )
    buf.extend(padding)
bookkeeping.misc += len(padding)
buf.tofile(fout)
write_padlen(padding, fout)
write_magic(fout)
fout.flush()
#bookkeeping.show()
