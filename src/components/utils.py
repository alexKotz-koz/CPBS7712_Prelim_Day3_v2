import os
import json


def toKmers(k, sequence):
    kmers = [sequence[i : i + k] for i in range(len(sequence) - k + 1)]
    return kmers


# useage: dictKeys(<python dict>, <filename.extension>, <fileLocation> with trailing '/' (optional))
def dictKeys(dict: dict, filename: str, fileLocation=None):
    """
    Args:
        <python dict>
        <filename.extension>
        <fileLocation> with trailing '/' (optional)
    """
    cwd = os.getcwd()
    if fileLocation:
        os.path.join(cwd, fileLocation)
    else:
        fileLocation = cwd
    with open(os.path.join(fileLocation, filename), "w") as file:
        json.dump(dict, file)
