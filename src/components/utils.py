def toKmers(k, sequence):
    kmers = [sequence[i : i + k] for i in range(len(sequence) - k + 1)]
    return kmers