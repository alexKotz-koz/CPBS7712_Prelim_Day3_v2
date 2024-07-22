import statistics
import json
import os
import pandas as pd
import logging


class QualityControl:
    def __init__(self, biosample):
        self.biosample = biosample
        self.logDataDir = "./data/logs"
        os.makedirs(self.logDataDir, exist_ok=True)
        self.outputDataDir = "./data/output_data"
        os.makedirs(self.outputDataDir, exist_ok=True)

    # Phred and P Error calculations: https://www.drive5.com/usearch/manual/quality_score.html
    # Input: character to convert
    # Output: Q-score
    def asciiToPhred(self, character):
        return ord(character) - 33

    # Input: Q-score
    # Output: probability error
    def calculatePError(self, q):
        return round(10 ** (-q / 10.0), 5)

    # Input: biosample dictionary
    # Output: cleaned biosample dictionary, dataframe, and quality control information used in reports and logging
    def qualityControl(self):
        biosample = self.biosample
        cleanedBiosample = {}
        qualityControlReport = {}
        qscores = []
        perrorscores = []
        sequenceLengths = []
        lengthOriginalBiosample = len(biosample)
        for id, read in biosample.items():
            id = id
            sequence = read["sequence"]
            quality = read["quality"]

            qscores = [self.asciiToPhred(q) for q in quality]

            trimmedScores = []
            trimmedSequence = []

            half_index = len(qscores) // 2
            for index in range(half_index, len(qscores)):
                if qscores[index] >= 20:
                    trimmedScores.append(qscores[index])
                    trimmedSequence.append(sequence[index])

            qscores = qscores[: len(qscores) // 2] + trimmedScores
            sequence = sequence[: len(sequence) // 2] + "".join(trimmedSequence)

            sequenceLengths.append(len(sequence))

            readMedianQ = statistics.median(qscores)
            readMeanQ = statistics.mean(qscores)
            if readMedianQ >= 20:
                sequence = "".join(sequence)
                cleanedBiosample[id] = {"sequence": sequence}
                qscoreLen = len(qscores)
                qscores = []

            qualityControlReport[id] = {
                "meanQ": readMeanQ,
                "medianQ": readMedianQ,
                "length": len(sequence),
            }
        avgSeqLength = sum(sequenceLengths) / len(sequenceLengths)

        qcReportFile = os.path.join(self.logDataDir, "QualityControlReport.json")
        with open(qcReportFile, "w") as file:
            json.dump(qualityControlReport, file)

        cleanedBiosampleFile = os.path.join(self.outputDataDir, "CleanedBioSample.json")
        with open(cleanedBiosampleFile, "w") as file:
            json.dump(cleanedBiosample, file)

        biosampleList = [{"id": id, **data} for id, data in cleanedBiosample.items()]
        biosampleDf = pd.DataFrame(biosampleList)

        logging.info("Quality Control: ")
        logging.info(f"\t# of cleaned reads: {len(cleanedBiosample)}")
        logging.info(f"\tAverage read length: {avgSeqLength}")
        logging.info(f"\tMinimum read length: {min(sequenceLengths)}")
        logging.info(f"\tMaximum read length: {max(sequenceLengths)}")
        qcMetaData = {
            "lengthOriginalBiosample": lengthOriginalBiosample,
            "lengthCleanedBiosample": len(cleanedBiosample),
            "averageReadLength": avgSeqLength,
            "minimumReadLength": min(sequenceLengths),
            "maximumReadLength": max(sequenceLengths),
        }

        return biosampleDf, min(sequenceLengths), qualityControlReport, qcMetaData
