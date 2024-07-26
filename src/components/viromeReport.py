import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import json
import math
import warnings
from components import utils

# ignoring deprecated feature in seaborn
warnings.filterwarnings("ignore", category=FutureWarning)


class ViromeReport:
    def __init__(
        self,
        contigs,
        virusesInBiosample,
        biosampleFile,
        biosampleFileLocation,
        qcMetadata,
    ):
        self.contigs = contigs
        self.virusesInBiosample = virusesInBiosample  # contains the viruses tested against the biosample (i.e. 4 for bat) ... can access contigsInVirus
        self.biosampleFile = biosampleFile
        self.biosampleFileLocation = biosampleFileLocation
        self.qcMetadata = qcMetadata
        self.reportDir = "data/reports"
        scriptDir = os.path.dirname(os.path.dirname(__file__))
        self.dataDir = os.path.join(scriptDir, "data")
        self.taxDir = os.path.join(self.dataDir, "taxonomy")
        os.makedirs(self.reportDir, exist_ok=True)

    # Input: sequence file
    # Ouput: number of reads in input file
    # Util fx
    def getNumReads(self, fileLocation):
        with open(fileLocation, "r") as file:
            lines = file.readlines()
            reads = [line for line in lines if line[0] == "@"]
            numReads = len(reads)
            return numReads

    # Input: viruses, biosample contigs
    # Output: virusAbundance dictionary,(# of contigs that match parts of virus)/total # of contigs
    def virusAbundance(self):
        virusAbundance = {}
        contigs = self.contigs
        totalNumContigs = len(contigs)
        virusesInBiosample = self.virusesInBiosample
        for virus in virusesInBiosample:
            vName = virus["virus"]
            numContigs = virus["numContigsInVirus"]

            vProportion = numContigs / totalNumContigs
            virusAbundance[vName] = {
                "abundance": vProportion * 100,
                "numContigsInVirus": numContigs,
            }
        return virusAbundance

    # Input: biosample dataframe (containing qc data), biosample name, pdf file name (for output)
    # Output: pdf report on qc for this run
    def createPDF(self, biosampleDf, biosampleName, ouputFileName):
        with PdfPages(ouputFileName) as pdf:
            fig, ax = plt.subplots(figsize=(12, 4))
            ax.axis("tight")
            ax.axis("off")
            plt.title(f"Biosample Information: {biosampleName}", fontsize=14, pad=20)

            biosample_table = ax.table(
                cellText=biosampleDf.values, colLabels=biosampleDf.columns, loc="center"
            )
            biosample_table.auto_set_font_size(False)
            biosample_table.set_fontsize(10)
            biosample_table.scale(1.2, 1.2)
            pdf.savefig(fig, bbox_inches="tight")

    # Input: virus dataframe containing taxonomy data
    # Output: Addition of the relative abundance column to the dataframe
    # Util fx
    def addColumnsToDf(self, virusDf):
        virusDf.set_index("virusName", inplace=True)
        virusDf["Relative Viral Abundance"] = np.nan
        virusDf["# Biosample Contigs in Virus"] = np.nan
        for virus, abundance in self.virusAbundance().items():
            if virus in virusDf.index:
                virusDf.loc[virus, "# Biosample Contigs in Virus"] = abundance[
                    "numContigsInVirus"
                ]
                virusDf.loc[virus, "Relative Viral Abundance"] = abundance["abundance"]

    # Input: constructor variables
    # Output: report files
    def generateReport(self):
        reportFile = "virome_report.txt"
        biosampleFileName = self.biosampleFile.replace(".fastq", "")
        # QC metadata
        lengthOriginalBiosample = self.qcMetadata["lengthOriginalBiosample"]
        lengthCleanedBiosample = self.qcMetadata["lengthCleanedBiosample"]
        avgReadLength = self.qcMetadata["averageReadLength"]
        minReadLen = self.qcMetadata["minimumReadLength"]
        maxReadLen = self.qcMetadata["maximumReadLength"]

        # Tax table
        if "bat" in self.biosampleFile:
            taxFile = os.path.join(self.taxDir, "CoV-Taxonomy.json")
            batVirusDf = pd.read_json(taxFile)

            self.addColumnsToDf(virusDf=batVirusDf)

            batVirusDf = batVirusDf[batVirusDf["Relative Viral Abundance"].notna()]

            batVirusDf.rename(index={"virusName": "Virus Name"}, inplace=True)
            batVirusDf.to_csv(
                os.path.join(self.reportDir, f"{biosampleFileName}_tax.csv"),
                index_label="Virus Name",
            )
        elif "biofilm" in self.biosampleFile or "cryoconite" in self.biosampleFile:
            taxFile = os.path.join(self.taxDir, "NCLDVGenes-Taxonomy.json")
            NCLDVGeneDf = pd.read_json(taxFile)

            self.addColumnsToDf(virusDf=NCLDVGeneDf)

            NCLDVGeneDf = NCLDVGeneDf[NCLDVGeneDf["Relative Viral Abundance"].notna()]

            NCLDVGeneDf.rename(index={"virusName": "Virus Name"}, inplace=True)
            NCLDVGeneDf.to_csv(
                os.path.join(self.reportDir, f"{biosampleFileName}_tax.csv"),
                index_label="Virus Name",
            )
        else:
            taxFile = os.path.join(self.taxDir, "synthetic-Taxonomy.json")
            print(f"Tax File: {taxFile}")
            synthDf = pd.read_json(taxFile)

            self.addColumnsToDf(virusDf=synthDf)

            synthDf = synthDf[synthDf["Relative Viral Abundance"].notna()]

            synthDf.rename(index={"virusName": "Virus Name"}, inplace=True)
            synthDf.to_csv(
                os.path.join(self.reportDir, f"{biosampleFileName}_tax.csv"),
                index_label="Virus Name",
            )

        # Add the total number of reads from the original biosample file to biosampleInfo.
        # This is used in the collect_results.sh file to calculate the percentage of reads tested in the expirement
        with open(os.path.join(self.dataDir, "numReads.json"), "r") as file:
            biosampleNumReads = json.load(file)
        ogBat = list(biosampleNumReads.items())[0][1]
        ogBiofilm = list(biosampleNumReads.items())[1][1]
        ogCryoconite = list(biosampleNumReads.items())[2][1]
        ogSynthetic = list(biosampleNumReads.items())[3][1]

        if "bat" in biosampleFileName:
            numOgReads = ogBat
        elif "biofilm" in biosampleFileName:
            numOgReads = ogBiofilm
        elif "cryoconite" in biosampleFileName:
            numOgReads = ogCryoconite
        elif "synthetic" in biosampleFileName:
            numOgReads = ogSynthetic

        biosampleInfo = {
            "Metric": [
                "Number of Reads in Biosample File",
                "Number of Reads in Cleaned Biosample File",
                "Average Read Length after Quality Control (bp)",
                "Minimum Read Length after Quality Control (bp)",
                "Maximum Read Length after Quality Control (bp)",
                "Number of Reads in Original Biosample File",
            ],
            "Value": [
                lengthOriginalBiosample,
                lengthCleanedBiosample,
                avgReadLength,
                minReadLen,
                maxReadLen,
                numOgReads,
            ],
        }

        biosampleDf = pd.DataFrame(biosampleInfo)

        self.createPDF(
            biosampleDf,
            biosampleFileName,
            os.path.join(self.reportDir, f"{biosampleFileName}.pdf"),
        )

        # Create report file
        reportFileLocation = os.path.join(self.reportDir, reportFile)
        with open(reportFileLocation, "a") as file:
            file.write("Biosample Information:\n")
            for metric, value in zip(biosampleInfo["Metric"], biosampleInfo["Value"]):
                file.write(f"\t{metric}: {value}\n")
            file.write("\nVirus Abundance in the biosample:\n")
            json.dump(self.virusAbundance(), file, indent=4)
