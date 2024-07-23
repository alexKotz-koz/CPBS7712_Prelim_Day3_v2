import pandas as pd
import logging
from logging.handlers import RotatingFileHandler

import queue
import argparse
import time
import os
import sys
import cProfile
import json

from components.importBioSample import ImportBioSample
from components.importVirus import ImportVirus
from components.qc import QualityControl

from components.readsToKmers import ReadsToKmers
from components.deBruijnGraph import DeBruijnGraph
from components.createContigs import CreateContigs

# Uncomment these lines to test the implementation of the smith-waterman algorithms and adjust the class instance initalization on line 183
# from components.searchForViruses_SW import SearchForViruses
# from components.searchForViruses_SW_PP import SearchForViruses
from components.searchForViruses import SearchString

from components.viromeReport import ViromeReport
from components.codeReport import CodeReport

logDir = "data/logs"
os.makedirs(logDir, exist_ok=True)
# Set up logging with a rotating file handler for the main.py file
log_file = os.path.join(logDir, "app.log")

handler = RotatingFileHandler(
    log_file,
    maxBytes=1024 * 1024,
    backupCount=5,
)

formatter = logging.Formatter("%(message)s")
handler.setFormatter(formatter)

logger = logging.getLogger()
logger.setLevel(logging.INFO)
logger.addHandler(handler)


def main():
    logging.info("Main: ")
    scriptDir = os.path.dirname(__file__)
    dataDir = "data"
    dataDir = os.path.join(scriptDir, "data")

    # arg setup and management
    parser = argparse.ArgumentParser(
        description="Metagenomic-based virome characterizer"
    )
    parser.add_argument(
        "-biosample", type=str, help="Fastq biosample file", required=True
    )
    parser.add_argument("-k", type=int, help="size of kmer", required=True)

    args = parser.parse_args()

    biosampleFile = args.biosample
    k = args.k
    logging.info(f"\tBioSample File: {biosampleFile}")
    logging.info(f"\tSize of K = {k}")

    componentRunTimes = {}

    viruses = {}
    virusDataDir = "./data/virus_data"
    RmYN02PrimersDataDir = "./data/virus_data/RmYN02Primers"
    

    # Glacial
    virusFile = "sequences_20240607_3345067.fasta"
    virusFile2 = "sequences_20240607_570283.fasta"
    virusFile3 = "sequences_20240607_9774926.fasta"
    virusFile4 = "sequences_20240607_5959983.fasta"
    ncbi4VirusesFile = "ncbi_random_4_viruses.fasta"
    NCLDVFile = "sequences_ Nucleocytoviricota.fasta"
    PolBGeneFile = "kegg_polB.fasta"
    A32GeneFile = "ncbi_A32.fna"
    D5GeneFile = "ncbi_D5.fna"
    RNAplGeneFile = "kegg_RNApl.fasta"
    RNApsGeneFile = "kegg_RNAps.fasta"
    mRNAcGeneFile = "kegg_mRNAc.fasta"
    RNRSFIIGeneFile = "kegg_RNRSFII.fasta"
    VLTF3GeneFile = "kegg_VLTF3.fasta"

    # Bat
    sarsCoV2File = "SarsCoV2.fasta"
    RmYN02FPrimer1 = "RmYN02FPrimer1.fasta"
    RmYN02RPrimer1 = "RmYN02RPrimer1.fasta"
    RmYN02FPrimer2 = "RmYN02FPrimer2.fasta"
    RmYN02RPrimer2 = "RmYN02RPrimer2.fasta"
    RmYN02FPrimer3 = "RmYN02FPrimer3.fasta"
    RmYN02RPrimer3 = "RmYN02RPrimer3.fasta"
    RmYN02FPrimer4 = "RmYN02FPrimer4.fasta"
    RmYN02RPrimer4 = "RmYN02RPrimer4.fasta"
    RmYN02FPrimer5 = "RmYN02FPrimer5.fasta"
    RmYN02RPrimer5 = "RmYN02RPrimer5.fasta"
    RmYN02FPrimer6 = "RmYN02FPrimer6.fasta"
    RmYN02RPrimer6 = "RmYN02RPrimer6.fasta"
    RmYN02FPrimer7 = "RmYN02FPrimer7.fasta"
    RmYN02RPrimer7 = "RmYN02RPrimer7.fasta"
    RmYN02FPrimer8 = "RmYN02FPrimer8.fasta"
    RmYN02RPrimer8 = "RmYN02RPrimer8.fasta"
    RmYN02Probe = "RmYN02Probe.fasta"

    

    # list of file locations
    NCLDVFileLocation = [os.path.join(virusDataDir, NCLDVFile)]
    ncbi4VirusesFileLocation = [os.path.join(virusDataDir, ncbi4VirusesFile)]
    NCLDVGeneFileLocations = [
        os.path.join(virusDataDir, PolBGeneFile),
        os.path.join(virusDataDir, A32GeneFile),
        os.path.join(virusDataDir, D5GeneFile),
        os.path.join(virusDataDir, RNAplGeneFile),
        os.path.join(virusDataDir, RNApsGeneFile),
        os.path.join(virusDataDir, mRNAcGeneFile),
        os.path.join(virusDataDir, RNRSFIIGeneFile),
        os.path.join(virusDataDir, VLTF3GeneFile),
    ]
    allGlacialVirusDataFileLocations = [
        os.path.join(virusDataDir, virusFile),
        os.path.join(virusDataDir, virusFile2),
        os.path.join(virusDataDir, virusFile3),
        os.path.join(virusDataDir, virusFile4),
        os.path.join(virusDataDir, NCLDVFile),
    ]
    batVirusFileLocations = [
        #os.path.join(virusDataDir, sarsCoV2File),
        os.path.join(RmYN02PrimersDataDir, RmYN02FPrimer1),
        os.path.join(RmYN02PrimersDataDir, RmYN02RPrimer1),
        os.path.join(RmYN02PrimersDataDir, RmYN02FPrimer2),
        os.path.join(RmYN02PrimersDataDir, RmYN02RPrimer2),
        os.path.join(RmYN02PrimersDataDir, RmYN02FPrimer3),
        os.path.join(RmYN02PrimersDataDir, RmYN02RPrimer3),
        os.path.join(RmYN02PrimersDataDir, RmYN02FPrimer4),
        os.path.join(RmYN02PrimersDataDir, RmYN02RPrimer4),
        os.path.join(RmYN02PrimersDataDir, RmYN02FPrimer5),
        os.path.join(RmYN02PrimersDataDir, RmYN02RPrimer5),
        os.path.join(RmYN02PrimersDataDir, RmYN02FPrimer6),
        os.path.join(RmYN02PrimersDataDir, RmYN02RPrimer6),
        os.path.join(RmYN02PrimersDataDir, RmYN02FPrimer7),
        os.path.join(RmYN02PrimersDataDir, RmYN02RPrimer7),
        os.path.join(RmYN02PrimersDataDir, RmYN02FPrimer8),
        os.path.join(RmYN02PrimersDataDir, RmYN02RPrimer8),
        os.path.join(RmYN02PrimersDataDir, RmYN02Probe)
    ]
    syntheticVirusFileLocation = [
        os.path.join(os.path.join(dataDir, "synthetic_data"), "synthetic_virus.fasta")
    ]

    if "synthetic" in biosampleFile:
        virusDataFileLocations = syntheticVirusFileLocation
    else:
        virusDataFileLocations = batVirusFileLocations

    if biosampleFile == "synthetic":
        biosampleFile = "synthetic_biosample.fastq"

    #### Main ####
    bioStart = time.time()
    importBioSampleInstance = ImportBioSample(biosampleFile=biosampleFile)
    biosample = importBioSampleInstance.importBioSample()
    bioStop = time.time()
    bioTotal = bioStop - bioStart
    componentRunTimes["importBioSample"] = bioTotal

    virusStart = time.time()
    importVirusInstance = ImportVirus()
    viruses = importVirusInstance.importVirusData(fileLocations=virusDataFileLocations)
    virusStop = time.time()
    virusTotal = virusStop - virusStart
    componentRunTimes["importVirus"] = virusTotal

    qcStart = time.time()
    qualityControlInstance = QualityControl(biosample=biosample)
    cleanedBiosample, minimumReadLength, qualityControlReport, qcMetadata = (
        qualityControlInstance.qualityControl()
    )
    qcStop = time.time()
    qcTotal = qcStop - qcStart
    componentRunTimes["qc"] = qcTotal

    if k > (minimumReadLength - 2):
        print(
            f"\nK must be at least one less than the size of the smallest read.\nMinimum read length for this sample is: {minimumReadLength}"
        )
        sys.exit(1)

    rtkStart = time.time()
    readsToKmersInstance = ReadsToKmers(readsData=cleanedBiosample, k=k)
    kmerPool = readsToKmersInstance.extractKmers()
    rtkStop = time.time()
    rtkTotal = rtkStop - rtkStart
    componentRunTimes["readsToKmers"] = rtkTotal
    logging.info(f"Time Stamp: Reads to Kmers finished in {rtkTotal}")
    print(f"Time Stamp: Reads to Kmers finished in {rtkTotal}")

    # Do not delete - used in search for viruses, added to managing logging with parallel processing
    with open("data/logs/r-kmerPool.json", "w") as file:
        json.dump(kmerPool, file)

    dbgStart = time.time()
    debruijnGraphInstance = DeBruijnGraph(kmerPool=kmerPool, k=k)
    nodes, edges = debruijnGraphInstance.constructGraph()
    dbgStop = time.time()
    dbgTotal = dbgStop - dbgStart
    componentRunTimes["deBruijnGraph"] = dbgTotal
    logging.info(f"Time Stamp: DeBruijn Graph finished in {dbgTotal}")
    print(f"Time Stamp: DeBruijn Graph finished in {dbgTotal}")

    ccStart = time.time()
    createContigsInstance = CreateContigs(graph=edges)
    contigs = createContigsInstance.createContigs()
    ccStop = time.time()
    ccTotal = ccStop - ccStart
    componentRunTimes["createContigs"] = ccTotal
    logging.info(f"Time Stamp: Create Contigs finished in {ccTotal}")
    print(f"Time Stamp: Create Contigs finished in {ccTotal}")

    # Modify the class instance call for testing the smith-waterman implementations
    sfvStart = time.time()
    searchForVirusesInstance = SearchString(
        viruses, "data/logs/r-kmerPool.json", contigs, k
    )
    virusesInBiosample = searchForVirusesInstance.searchString()
    sfvStop = time.time()
    sfvTotal = sfvStop - sfvStart
    componentRunTimes["searchForViruses"] = sfvTotal
    logging.info(f"Time Stamp: Find Viruses finished in {sfvTotal}")
    print(f"Time Stamp: Find Viruses finished in {sfvTotal}")

    #viromeReportInstance = ViromeReport(
    #    contigs, virusesInBiosample, biosampleFile, qcMetadata
    #)
    #viromeReportInstance.generateReport()


if __name__ == "__main__":
    start = time.time()
    main()
    stop = time.time()
    logging.info(f"Project execution time: {stop-start}")
