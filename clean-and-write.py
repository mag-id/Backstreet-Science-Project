
import sys
import projectlib as my

'''clean-and-write script'''

# directories
template, folder = str(sys.argv[1]), str(sys.argv[2])

# initialization
pattern = r'\s+([a-zA-Z]*)\s+(\-?\d+\.\d+)\s+\d\s+(\-?\d+\.\d+)\s+\d\s+(\-?\d+\.\d+)\s+\d\s+'
columns = 'atom_name', 'x', 'y', 'z'
template_xyz = my.xyzextract(template, pattern=pattern, columns=columns)

for path in my.find(folder, pattern=r'.+\.mop'):

    xyz = my.xyzextract(path, pattern=pattern, columns=columns)
    # Check and drop
    for string in my.atomdrop(template_xyz, xyz):
        xyz = xyz.drop(string)
    # Rewrite
    my.xyzrewrite(xyz, path, pattern=pattern)
