import os
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import warnings

# ignoring deprecated feature in seaborn
warnings.filterwarnings("ignore", category=FutureWarning)


class CodeReport:
    def __init__(
        self,
        qualityControlReport,
        contigs,
        componentRunTimes,
        qcMetadata,
        biosampleFile,
        k,
    ):
        self.qcReport = qualityControlReport
        self.contigs = contigs
        self.componentRunTimes = componentRunTimes
        self.qcMetadata = qcMetadata
        self.biosampleFile = biosampleFile
        self.k = k
        self.reportDir = "data/reports"
        os.makedirs(self.reportDir, exist_ok=True)

    # Contigs generated
    # Size of k
    # Time for each component completion

    def qcReportAnalysis(self):
        qcReportDf = pd.DataFrame(self.qcReport).T

        fig, axs = plt.subplots(3, figsize=(10, 10))

        yLabels = [
            "Mean Q-Score (per read)",
            "Median Q-Score (per read)",
            "Length (per read)",
        ]
        colors = ["blue", "green", "red"]
        for index, column in enumerate(qcReportDf.columns):
            sns.lineplot(
                data=qcReportDf[column], ax=axs[index], marker="o", color=colors[index]
            )
            axs[index].set_xticklabels([])
            axs[index].set_ylabel(yLabels[index], weight="bold")

        fig.text(0.5, 0.04, "Reads", ha="center", va="center", weight="bold")

        fig.suptitle("Quality Control Report", weight="bold")
        plt.subplots_adjust(bottom=0.1)

        figFile = os.path.join(self.reportDir, "QualityControl.png")
        plt.savefig(figFile)

    def assemblyStatistics(self):
        # https://www.ncbi.nlm.nih.gov/pmc/articles/PMC11020464/
        # https://www.nature.com/articles/s41467-022-35735-y fig 1

        contigs = self.contigs
        numContigs = len(contigs)
        contigLens = []
        for contig in self.contigs:
            contigLens.append(len(contig))
        minContig = min(contigLens)
        maxContig = max(contigLens)
        avgContig = sum(contigLens) / len(contigLens)

        return numContigs, minContig, maxContig, avgContig

    def generateReport(self):
        reportFile = "code_report.txt"

        numReads = self.qcReportAnalysis()
        numContigs, minContig, maxContig, avgContig = self.assemblyStatistics()

        reportFileLocation = os.path.join(self.reportDir, reportFile)
        with open(reportFileLocation, "w") as file:
            file.write("Code Report: \n\n")
            file.write(f"Biosample File: {self.biosampleFile}\n")
            file.write(f"Size of User-defined K: {self.k}\n\n")
            file.write("Component Execution Times: \n")
            for name, time in self.componentRunTimes.items():
                file.write(f"\t{name}.py completed in {time} seconds\n")
            file.write("\nContig Information: \n")
            file.write(f"\tTotal Number of contigs created: {numContigs}\n")
            file.write(f"\tMinimum Contig Length: {minContig} bp's\n")
            file.write(f"\tMaximum Contig Length: {maxContig} bp's\n")
            file.write(f"\tAverage Contig Length: {avgContig} bp\n")
