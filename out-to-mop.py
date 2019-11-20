
import sys
import projectlib as my

'''out-to-mop script'''

folder = str(sys.argv[1])

mop_path = my.construct(folder, '/SP/_MOPs')
csv_path = my.construct(folder, '/SP/_CSVs')
my.create(mop_path, csv_path)

for path in my.find(folder, pattern=r'test\d+\.out'):

    stem, suffix = my.deconstruct(path)[1, 2]
    comment = f'Conformer {stem}.{suffix}'

    my.xyz2mop(my.xyzextract(path), my.charged(path), stem, comment, mop_path, csv_path)
