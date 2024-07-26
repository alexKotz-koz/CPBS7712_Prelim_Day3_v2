import pandas as pd
from collections import defaultdict
from components import utils


class ReadsToKmers:
    def __init__(self, readsData, k):
        self.readsData = readsData
        self.k = k

    # Input: reads from sample
    # Output: k-mers from reads. Each k-mer has an id.
    def extractKmers(self):
        kmerPool = defaultdict(lambda: defaultdict(list))
        readsData = self.readsData
        k = self.k

        for read in readsData.itertuples():
            id = read.id
            sequence = read.sequence
            kmers = utils.toKmers(k, sequence)

            for index, kmer in enumerate(kmers):
                kmerPool[kmer][id].append({index: index + k})

        return kmerPool
