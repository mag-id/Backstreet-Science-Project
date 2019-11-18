
import pathlib
import json


def deconstruct(path: str) -> (str, str, str):
    """"""
    pop = pathlib.Path(path)

    folder = pop.parent
    stem = pop.stem
    suffix = pop.suffix

    return str(folder), stem, suffix


def construct(folder: str = pathlib.Path.cwd(), stem: str = '', suffix: str = '') -> str:
    """"""
    folder = pathlib.Path(folder)
    name = pathlib.Path(stem + suffix)

    return str(pathlib.Path.joinpath(folder, name))


def create(path: str) -> True:
    """"""
    path = pathlib.Path(path)
    pathlib.Path.mkdir(path, parents=True, exist_ok=False)
    return True


def json_save(path: str, dictionary: dict) -> True:
    """"""
    with open(path, 'a') as f:
        json.dump(dictionary, f)
    return True


def json_load(path: str) -> dict:
    """"""
    with open(path) as f:
        dictionary = json.load(f)
    return dictionary
