bases_and_eof = ['a', 'c', 'g', 't', 'n', 'eof']
base_codes = dict(zip(bases_and_eof, range(6)))

# truncated binary encoding of single base
bases_tbe_codes = ['00', '01', '100', '101', '110', '111']
base2tbe = dict(zip(bases_and_eof, bases_tbe_codes))
tbe2base = dict(zip(bases_tbe_codes, bases_and_eof))
