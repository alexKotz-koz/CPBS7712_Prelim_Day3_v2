import pandas as pd
from collections import defaultdict


class ReadsToKmers:
    def __init__(self, readsData, k):
        self.readsData = readsData
        self.k = k

    def extractKmers(self):
        kmerPool = defaultdict(lambda: defaultdict(list))
        readsData = self.readsData
        k = self.k

        for read in readsData.itertuples():
            id = read.id
            sequence = read.sequence
            kmers = [sequence[i : i + k] for i in range(len(sequence) - k + 1)]

            for index, kmer in enumerate(kmers):
                kmerPool[kmer][id].append({index: index + k})

        return kmerPool
