import os
import sys
import json
import pandas as pd
import numpy as np
import logging
import time


class SearchForViruses:
    def __init__(self, viruses, contigs, k):
        self.viruses = viruses
        self.contigs = contigs
        self.k = k
        self.logDataDir = "./data/logs"
        self.outputDataDir = "./data/output_data"

    def smith_waterman(
        self, contig, virus, matchScore=3, mismatchScore=-1, gapPenalty=-2
    ):
        # initialize scoring matrix
        rows = len(contig) + 1
        cols = len(virus) + 1
        scoreMatrix = np.zeros((rows, cols), dtype=int)

        # initialize traceback matrix
        tracebackMatrix = np.zeros((rows, cols), dtype=int)

        # fill in scoring and traceback matrices
        for i in range(1, rows):
            for j in range(1, cols):
                match = scoreMatrix[i - 1][j - 1] + (
                    matchScore if contig[i - 1] == virus[j - 1] else mismatchScore
                )
                delete = scoreMatrix[i - 1][j] + gapPenalty
                insert = scoreMatrix[i][j - 1] + gapPenalty
                scoreMatrix[i][j] = max(match, delete, insert, 0)
                tracebackMatrix[i][j] = np.argmax([0, match, delete, insert])

        # find the maximum score in the matrix
        maxScore = np.max(scoreMatrix)

        # find the index of the maximum score
        maxIndex = np.unravel_index(np.argmax(scoreMatrix), scoreMatrix.shape)

        # traceback to reconstruct the aligned sequences
        alignedContig = ""
        alignedVirus = ""
        i, j = maxIndex
        startPosition = j
        while i > 0 and j > 0 and scoreMatrix[i][j] != 0:
            if tracebackMatrix[i][j] == 1:
                alignedContig = contig[i - 1] + alignedContig
                alignedVirus = virus[j - 1] + alignedVirus
                i -= 1
                j -= 1
            elif tracebackMatrix[i][j] == 2:
                alignedContig = contig[i - 1] + alignedContig
                alignedVirus = "-" + alignedVirus
                i -= 1
            else:
                alignedContig = "-" + alignedContig
                alignedVirus = virus[j - 1] + alignedVirus
                j -= 1
        endPosition = j
        return maxScore, alignedContig, alignedVirus, startPosition, endPosition

    def calculateCoverage(self, alignments, virusLength):
        coverageArray = np.zeros(virusLength, dtype=bool)
        for alignment in alignments:
            start = alignment["startPosition"]
            end = alignment["endPosition"]
            if start != -1 and end != -1:
                coverageArray[start:end] = True
        coverage = np.sum(coverageArray) / virusLength * 100
        return coverage

    def search(self):
        logging.info("Search for Viruses: ")
        virusContigObj = {}
        viruses = self.viruses
        vIndex = 0
        for virus, virusData in viruses.items():
            vStart = time.time()
            virusSequence = virusData["sequence"]
            vIndex += 1
            alignments = []
            for index, contig in enumerate(self.contigs):
                maxScore, alignedContig, alignedVirus, startPosition, endPosition = (
                    self.smith_waterman(contig, virusSequence)
                )
                print(
                    f"Length of aligned contig: {len(alignedContig)}. Length of aligned Virus:{len(alignedVirus)}"
                )
                if virusData["name"] not in virusContigObj:
                    virusContigObj[virusData["name"]] = [
                        {
                            "contig": index + 1,
                            "alignmentScore": int(maxScore),
                            "alignedContigSubstring": alignedContig,
                            "alignedVirusSubstring": alignedVirus,
                            "lengthOfAlignment": int(len(alignedContig)),
                            "startPosition": int(startPosition),
                            "endPosition": int(endPosition),
                        }
                    ]
                else:
                    virusContigObj[virusData["name"]].append(
                        {
                            "contig": index + 1,
                            "alignmentScore": int(maxScore),
                            "alignedContigSubstring": alignedContig,
                            "alignedVirusSubstring": alignedVirus,
                            "lengthOfAlignment": int(len(alignedContig)),
                            "startPosition": int(startPosition),
                            "endPosition": int(endPosition),
                        }
                    )

            vStop = time.time()
            logging.info(f"\tVirus {virusData['name']} completed in {vStop-vStart}")
            coverage = self.calculateCoverage(alignments, len(virusSequence))
            virusContigObj[virusData["name"]].append({"coverage": coverage})
            logging.info(f"\tVirus coverage: {coverage}")
        virusContigObjFile = os.path.join(self.outputDataDir, "VirusContigObject.json")
        with open(virusContigObjFile, "w") as file:
            json.dump(virusContigObj, virusContigObjFile)
