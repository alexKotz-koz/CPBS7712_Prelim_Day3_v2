import os
import sys
import json
import pandas as pd
import numpy as np
import time

from concurrent.futures import ProcessPoolExecutor


class SearchForViruses:
    def __init__(self, viruses, contigs, k):
        self.viruses = viruses
        self.contigs = contigs
        self.k = k
        self.outputDataDir = "./data/output_data"

    def smith_waterman(
        self, contig, virus, matchScore=3, mismatchScore=-1, gapPenalty=-2
    ):
        # initialize scoring matrix
        rows = len(contig) + 1
        cols = len(virus) + 1
        scoreMatrix = np.zeros((rows, cols), dtype=int)

        # initialize traceback matrix
        traceBackMatrix = np.zeros((rows, cols), dtype=int)

        # starting from the second row and second column,
        for i in range(1, rows):
            for j in range(1, cols):
                match = scoreMatrix[i - 1][j - 1] + (
                    matchScore if contig[i - 1] == virus[j - 1] else mismatchScore
                )
                delete = scoreMatrix[i - 1][j] + gapPenalty
                insert = scoreMatrix[i][j - 1] + gapPenalty
                scoreMatrix[i][j] = max(match, delete, insert, 0)
                traceBackMatrix[i][j] = np.argmax([0, match, delete, insert])

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
            if traceBackMatrix[i][j] == 1:
                alignedContig = contig[i - 1] + alignedContig
                alignedVirus = virus[j - 1] + alignedVirus
                i -= 1
                j -= 1
            elif traceBackMatrix[i][j] == 2:
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
        # create a bool array to store coverage
        coverageArray = np.zeros(virusLength, dtype=bool)
        for alignment in alignments:
            start = alignment["startPosition"]
            end = alignment["endPosition"]
            if start != -1 and end != -1:
                coverageArray[start:end] = True
        coverage = np.sum(coverageArray) / virusLength * 100
        return coverage

    def processContig(self, args):
        # generator function for processing each contig
        index, contig, virusSequence, virusData = args
        maxScore, alignedContig, alignedVirus, startPosition, endPosition = (
            self.smith_waterman(contig, virusSequence)
        )
        result = {
            "contig": index + 1,
            "alignmentScore": int(maxScore),
            "alignedContigSubstring": alignedContig,
            "alignedVirusSubstring": alignedVirus,
            "lengthOfAlignment": int(len(alignedContig)),
            "startPosition": int(startPosition),
            "endPosition": int(endPosition),
        }
        return virusData["name"], result

    def search(self):
        virusContigObj = {}
        viruses = self.viruses
        vIndex = 0

        totalContigs = len(self.contigs)
        processedContigs = 0

        for virus, virusData in viruses.items():
            print(f"Virus to align: {virus}")
            print(f"Length of virus: {len(virusData['sequence'])}")

            virusSequence = virusData["sequence"]
            vIndex += 1
            alignments = []

            with ProcessPoolExecutor() as executor:
                results = executor.map(
                    self.processContig,
                    [
                        (index, contig, virusSequence, virusData)
                        for index, contig in enumerate(self.contigs)
                    ],
                )
                for virusName, result in results:
                    processedContigs += 1
                    # print remaining 'time' for each virus
                    if processedContigs == totalContigs // 4:
                        print(
                            f"25% of the contigs have been processed in {round(time.time() - vStart,3)} seconds, from the start."
                        )
                    elif processedContigs == totalContigs // 2:
                        print(
                            f"50% of the contigs have been processed in {round(time.time() - vStart,3)} seconds, from the start."
                        )
                    elif processedContigs == (totalContigs * 3) // 4:
                        print(
                            f"75% of the contigs have been processed in {round(time.time() - vStart,3)} seconds, from the start."
                        )

                    if virusName not in virusContigObj:
                        virusContigObj[virusName] = [result]
                    else:
                        virusContigObj[virusName].append(result)
            coverage = self.calculateCoverage(alignments, len(virusSequence))
            virusContigObj[virusData["name"]].append({"coverage": coverage})
            print(f"\tVirus coverage: {coverage}")

        virusContigObjFile = os.path.join(self.outputDataDir, "VirusContigObject.json")
        with open(virusContigObjFile, "w") as file:
            json.dump(virusContigObj, file)
