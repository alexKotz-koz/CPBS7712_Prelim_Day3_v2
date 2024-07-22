import json
import os
import logging
import time

logDir = "data/logs"
os.makedirs(logDir, exist_ok=True)

# Set up logging for the searchForViruses.py file
search_logger = logging.getLogger(__name__)
search_logger.setLevel(logging.INFO)

search_handler = logging.FileHandler(os.path.join(logDir, "searchForViruses.log"))
formatter = logging.Formatter("%(asctime)s - %(name)s - %(levelname)s - %(message)s")
search_handler.setFormatter(formatter)

search_logger.addHandler(search_handler)


class SearchString:
    def __init__(self, viruses, readsKmerPoolFile, contigs, k):
        self.viruses = viruses
        self.contigs = contigs
        with open(readsKmerPoolFile, "r") as file:
            readsKmerPool = json.load(file)
        self.readsKmerPool = readsKmerPool
        self.k = k
        self.maxHammingDistance = 2

    # Input: virus sequence
    # Output: virus kmer pool
    def virusToKmers(self, sequence):
        kmerPool = {}
        for index in range(len(sequence) - self.k + 1):
            kmer = sequence[index : index + self.k]
            if kmer not in kmerPool:
                kmerPool[kmer] = [index]
            else:
                kmerPool[kmer].append(index)
        return kmerPool

    # Input: virus kmerpool
    def kmerPoolsToFile(self, virusKmerPool):
        with open("data/logs/v-kmerPool.json", "w") as file:
            json.dump(virusKmerPool, file)

    # Input: virus k-mer and contig k-mer
    # Output: hamming distance between the two k-mers
    def hammingDistance(self, virus, contig):
        distance = 0
        for index in range(len(virus)):
            if contig[index] != virus[index]:
                distance += 1
        return distance

    # Input: virus k-mer pool and virus object
    # Output: Contigs info object contain information about each contig k-mer that aligned (or didn't align) with each virus
    def createContigsInfo(self, virus, virusKmerPool):
        contigsInfo = []
        for id, contig in enumerate(self.contigs):
            contigLen = len(contig)
            contigKmers = {
                contig[i : i + self.k] for i in range(contigLen - self.k + 1)
            }
            kmerCount = 0
            contigInfo = {
                "contigId": id + 1,
                "contig": contig,
                "length": contigLen,
                "v-kmers": [],
            }
            for virusKmer in virusKmerPool:
                for contigKmer in contigKmers:
                    distance = self.hammingDistance(virusKmer, contigKmer)
                    if distance <= self.maxHammingDistance:
                        kmerCount += 1
                        contigInfo["v-kmers"].append(
                            {
                                virusKmer: {
                                    "indexOfVKmerInContig": [
                                        contig.index(contigKmer),
                                        contig.index(contigKmer) + self.k,
                                    ],
                                    "hammingDistance": distance,
                                    "kmerLength": len(contigKmer),
                                }
                            }
                        )
            contigInfo["kmerCount"] = kmerCount
            contigsInfo.append(contigInfo)
        return contigsInfo

    # Input: virus and contig objects
    # Output: information about the contigs that aligned (or didn't align) with each virus
    def searchString(self):
        logging.info("Search For Viruses:\n")
        virusesInBiosample = []
        # break all viruses into kmers
        virusKmerPools = {
            virus_id: self.virusToKmers(virus["sequence"])
            for virus_id, virus in self.viruses.items()
        }

        for virusId, virus in self.viruses.items():
            vStart = time.time()
            contigsExistInVirus = []
            contigsTested = []
            virusKmerPool = virusKmerPools[virusId]

            contigsInfo = self.createContigsInfo(virus["name"], virusKmerPool)

            for contig in contigsInfo:
                if contig["kmerCount"] > 0:
                    contigsExistInVirus.append(contig)
                contigsTested.append(contig)
            virusesInBiosample.append(
                {
                    "virus": virus["name"],
                    "numContigsInVirus": len(contigsExistInVirus),
                    "contigsInVirus": contigsExistInVirus,
                    "contigsTested": contigsTested,
                }
            )
            vStop = time.time()

            print(f"\t{len(contigsExistInVirus)} contigs align with {virus['name']} ")

            logging.info(
                f"\t{len(contigsExistInVirus)} contigs align with {virus['name']} "
            )
            logging.info(f"\tVirus {virus['name']} split time: {vStop-vStart}")
            logging.info(f"\tVirus: {virus['name']} length: {len(virus['sequence'])}bp")

        with open("data/output_data/virusesInBiosample.json", "w") as file:
            json.dump(virusesInBiosample, file)

        return virusesInBiosample
