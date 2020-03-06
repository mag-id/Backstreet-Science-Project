
from alt_gaocs import class_structure as my

test_struct = my.Structure.read('small_test.sdf')
print(test_struct.name, '\n', test_struct.coordinates, '\n', test_struct.graph)

test_struct.write('/Users/gamag/PycharmProjects/BSP/alt_gaocs/test_struct.mop')
