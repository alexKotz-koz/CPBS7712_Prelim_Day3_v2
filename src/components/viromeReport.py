import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import json
import math
import warnings

# ignoring deprecated feature in seaborn
warnings.filterwarnings("ignore", category=FutureWarning)


class ViromeReport:
    def __init__(
        self,
        contigs,
        virusesInBiosample,
        biosampleFile,
        qcMetadata,
    ):
        self.contigs = contigs
        self.virusesInBiosample = virusesInBiosample
        self.biosampleFile = biosampleFile
        self.qcMetadata = qcMetadata
        self.reportDir = "data/reports"
        scriptDir = os.path.dirname(os.path.dirname(__file__))
        dataDir = "data"
        dataDir = os.path.join(scriptDir, "data")
        self.taxDir = os.path.join(dataDir, "taxonomy")
        os.makedirs(self.reportDir, exist_ok=True)
        # print(f"Viruses in Biosample: {self.virusesInBiosample["contigsInVirus"]}")

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
            virusAbundance[vName] = {"abundance": vProportion * 100}
        return virusAbundance

    def createPDF(self, biosampleDf, biosampleName, ouputFileName):
        with PdfPages(ouputFileName) as pdf:
            fig, ax = plt.subplots(figsize=(12, 4))
            ax.axis("tight")
            ax.axis("off")
            plt.title(f"Biosample Information: {biosampleName}", fontsize=14, pad=20)

            # Create table for biosample information
            biosample_table = ax.table(
                cellText=biosampleDf.values, colLabels=biosampleDf.columns, loc="center"
            )
            biosample_table.auto_set_font_size(False)
            biosample_table.set_fontsize(10)
            biosample_table.scale(1.2, 1.2)
            pdf.savefig(fig, bbox_inches="tight")

    def addAbundaceToDf(self, virusDf):
        virusDf.set_index("virusName", inplace=True)
        virusDf["Relative Viral Abundance"] = np.nan

        # virusAubndance "{name, %contigs in virus}"
        for virus, abundance in self.virusAbundance().items():
            if virus in virusDf.index:
                virusDf.loc[virus, "Relative Viral Abundance"] = abundance["abundance"]

    def generateReport(self):
        reportFile = "virome_report.txt"
        biosampleFileName = self.biosampleFile.strip(".fastq")

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

            self.addAbundaceToDf(virusDf=batVirusDf)

            # Filter the DataFrame to include only the rows where 'Relative Viral Abundance' is not NaN
            batVirusDf = batVirusDf[batVirusDf["Relative Viral Abundance"].notna()]

            batVirusDf.rename(index={"virusName": "Virus Name"}, inplace=True)
            batVirusDf.to_csv(
                os.path.join(self.reportDir, f"{biosampleFileName}_tax.csv"),
                index_label="Virus Name",
            )
        elif "biofilm" in self.biosampleFile or "cryoconite" in self.biosampleFile:
            taxFile = os.path.join(self.taxDir, "NCLDVGenes-Taxonomy.json")
            NCLDVGeneDf = pd.read_json(taxFile)

            self.addAbundaceToDf(virusDf=NCLDVGeneDf)

            # Filter the DataFrame to include only the rows where 'Relative Viral Abundance' is not NaN
            NCLDVGeneDf = NCLDVGeneDf[NCLDVGeneDf["Relative Viral Abundance"].notna()]

            NCLDVGeneDf.rename(index={"virusName": "Virus Name"}, inplace=True)
            NCLDVGeneDf.to_csv(
                os.path.join(self.reportDir, f"{biosampleFileName}_tax.csv"),
                index_label="Virus Name",
            )
        else:
            taxFile = os.path.join(self.taxDir, "synthetic-Taxonomy.json")
            synthDf = pd.read_json(taxFile)

            self.addAbundaceToDf(virusDf=synthDf)

            # Filter the DataFrame to include only the rows where 'Relative Viral Abundance' is not NaN
            synthDf = synthDf[synthDf["Relative Viral Abundance"].notna()]

            synthDf.rename(index={"virusName": "Virus Name"}, inplace=True)
            synthDf.to_csv(
                os.path.join(self.reportDir, f"{biosampleFileName}_tax.csv"),
                index_label="Virus Name",
            )
        biosampleInfo = {
            "Metric": [
                "Number of Reads in Original Biosample File",
                "Number of Reads in Cleaned Biosample File",
                "Average Read Length after Quality Control (bp)",
                "Minimum Read Length after Quality Control (bp)",
                "Maximum Read Length after Quality Control (bp)",
            ],
            "Value": [
                lengthOriginalBiosample,
                lengthCleanedBiosample,
                avgReadLength,
                minReadLen,
                maxReadLen,
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
