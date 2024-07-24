import os
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import json
import math
import warnings

# ignoring deprecated feature in seaborn
warnings.filterwarnings("ignore", category=FutureWarning)


# number of contigs per virus with thier distance
class ViromeReport:
    def __init__(self, contigs, virusesInBiosample, biosampleFile, qcMetadata, virusDataFileLocations):
        self.contigs = contigs
        self.virusesInBiosample = virusesInBiosample
        self.biosampleFile = biosampleFile
        self.qcMetadata = qcMetadata
        self.virusDataFileLocations = virusDataFileLocations
        self.reportDir = "data/reports"
        scriptDir = os.path.dirname(os.path.dirname(__file__))
        dataDir = "data"
        dataDir = os.path.join(scriptDir, "data")
        self.taxDir = os.path.join(dataDir, "taxonomy")
        os.makedirs(self.reportDir, exist_ok=True)

    def virusAbundance(self):
        virusAbundance = {}
        contigs = self.contigs
        totalNumContigs = len(contigs)
        vInB = self.virusesInBiosample
        for virus in vInB:
            vName = virus["virus"]
            numContigs = virus["numContigsInVirus"]
            vProportion = numContigs / totalNumContigs
            virusAbundance[vName] = {"abundance": vProportion * 100}
        return virusAbundance

    def shannonDiversity(self, virusAbundance):
        # Shannon index diversity formula: https://www.statology.org/shannon-diversity-index/
        pITimeslnPis = []
        for index, virus in virusAbundance.items():
            if virus["abundance"] == 0:
                continue
            pI = virus["abundance"] / 100
            lnPi = math.log(pI)
            pITimeslnPi = pI * lnPi
            pITimeslnPis.append(pITimeslnPi)
        shannonDiversityIndex = sum(pITimeslnPis) * -1
        return shannonDiversityIndex

    def generateReport(self):
        reportFile = "virome_report.txt"

        virusAbundance = self.virusAbundance()
        df = pd.DataFrame(virusAbundance).T.reset_index()
        df.columns = ["Virus", "Abundance"]

        sns.barplot(
            x="Virus",
            y="Abundance",
            data=df,
        )
        plt.xticks(rotation=45, fontsize="x-small")
        plt.title(f"Virus Abundance in {self.biosampleFile}")
        figFile = os.path.join(self.reportDir, "VirusAbundance.png")
        plt.savefig(figFile, bbox_inches="tight")

        # shannon diversity index
        shannonIndex = self.shannonDiversity(virusAbundance)

        # qc metadata
        lengthOriginalBiosample = self.qcMetadata["lengthOriginalBiosample"]
        lengthCleanedBiosample = self.qcMetadata["lengthCleanedBiosample"]
        avgReadLength = self.qcMetadata["averageReadLength"]
        minReadLen = self.qcMetadata["minimumReadLength"]
        maxReadLen = self.qcMetadata["maximumReadLength"]

        #taxa table
        taxFile = os.path.join(self.taxDir, 'CoV-Taxonomy.json')
        sarsCov2TaxaDf = pd.read_json(taxFile)
        print(sarsCov2TaxaDf.head())
        # Set 'virusName' as the index
        sarsCov2TaxaDf.set_index('virusName', inplace=True)
        sarsCov2TaxaDf['abundance'] = sarsCov2TaxaDf.index.to_series().apply(
            lambda x: virusAbundance['Severe acute respiratory syndrome coronavirus 2']['abundance'] if 'SARS-CoV-2' in x 
            else (virusAbundance['Middle East respiratory syndrome-related coronavirus']['abundance'] if 'MERS-CoV' in x else '0')
        )


        sarsCov2TaxaDf.to_csv(os.path.join(self.reportDir, f"{self.biosampleFile}_tax.csv"))


        # create report file
        reportFileLocation = os.path.join(self.reportDir, reportFile)
        with open(reportFileLocation, "a") as file:
            file.write("\n\nVirome Report:\n\n")
            file.write("Biosample Information:\n")
            file.write(
                f"\tNumber of Reads in Original Biosample File: {lengthOriginalBiosample}\n"
            )
            file.write(
                f"\tNumber of Reads in Biosample after Quality Control: {lengthCleanedBiosample}\n\n"
            )
            file.write(
                f"\tAverage Read Length after Quality Control: {avgReadLength} bp\n"
            )
            file.write(
                f"\tMinimum Read Length after Quality Control: {minReadLen} bp\n"
            )
            file.write(
                f"\tMaximum Read Length after Quality Control: {maxReadLen} bp\n\n"
            )
            file.write(
                f"Shannon Diversity Index for {self.biosampleFile} community: {round(shannonIndex, 3)}\n\n"
            )

            file.write(f"Virus Abundance in {self.biosampleFile}:\n")
            json.dump(virusAbundance, file, indent=4)
