
import pybel
import filesys as fs


def confslice(inp_type: str, inp_path: str, out_type: str, out_path: str) -> True:
    """"""
    for number, structure in enumerate(pybel.readfile(inp_type, inp_path)):

        path, stem, suffix = fs.deconstruct(out_path)
        new_path = fs.construct(path, f'{stem}_{number}', suffix)

        conformer = pybel.Outputfile(out_type, new_path, overwrite=False)
        conformer.write(structure)

    return True
