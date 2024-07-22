import json
import logging
import time
from collections import defaultdict, deque
import os


class CreateContigs:
    def __init__(self, graph):
        self.graph = graph
        self.edgesCount = defaultdict(lambda: [0, 0])  # [incoming, outgoing]
        self.allPaths = []
        scriptDir = os.path.dirname(os.path.dirname(__file__))
        dataDir = "data"
        dataDir = os.path.join(scriptDir, "data")
        self.logsDataDir = os.path.join(dataDir, "logs")
        self.outputDataDir = os.path.join(dataDir, "output_data")
        os.makedirs(self.logsDataDir, exist_ok=True)
        os.makedirs(self.outputDataDir, exist_ok=True)

    # Input: edge list from readsToKmers
    # Output: a list of all start nodes (nodes that only have outgoing edges)
    def findStartNodes(self):
        startNodes = []

        for node, edges in self.graph.items():
            self.edgesCount[node][1] += len(edges)
            for edge in edges:
                self.edgesCount[edge][0] += 1

        startNodes = [
            node for node, counts in self.edgesCount.items() if counts[0] == 0
        ]
        return self.edgesCount, startNodes

    # Input: node to check
    # Output: bool, node is end of path or not
    def checkIfLastNode(self, currentNode):
        # if the currentNode is not in the graph or the # of outgoing edges for the current node is 0, return True
        return currentNode not in self.graph or self.edgesCount[currentNode][1] == 0

    # Input: node to check
    # Output: bool, list of children (nodes connected to current node via an outgoing edge from the current node)
    def lookForChildren(self, currentNode):
        if len(self.graph[currentNode]) > 1:
            return True, self.graph[currentNode]
        else:
            return False, []

    # Input: node that has children from original path, original path
    # Output: a continuation of the original path (can yeild more than one possible path)
    def followSubPath(self, startNode, originalPath):
        stack = deque([(startNode, originalPath + [startNode])])
        # continuation of modified DFS, creating a new path everytime a split occurs (2 or more outgoing edges from a single node)
        while stack:
            currentNode, path = stack.pop()
            if self.checkIfLastNode(currentNode=currentNode):
                yield path  # pause followSubPath execution, redirect control to followPath for consumption of yielded path
                continue

            children = self.graph.get(currentNode, [])

            for child in children:
                """
                - The uncommented code block below is not covering: loops (a path that starts and stops on the same node with n nodes in between), bubble (self repeating loop).
                - To handle one loop, replace the existing if statement with the commented if statement below.
                """
                # if path.count(child) < 2: #handles one loop
                if child not in path:  # collapses loop
                    newPath = path + [
                        child
                    ]  # create a new copy of path inside the loop
                    if self.checkIfLastNode(currentNode=child):
                        yield newPath  # pause followSubPath execution, redirect control to followPath for consumption of yielded newPath
                    else:
                        stack.append((child, newPath))

    # Input: start node and graph (edge list)
    # Output: allPaths object that contains all possible paths through the graph
    def followPath(self, startNode, visited=None):
        if visited is None:  # handle first iteration
            visited = deque([])
        unfinishedPaths = deque([])
        stack = deque([startNode])  # deque used to simulate a doubly linked list

        """
        This is a modified depth-first search (DFS). In this method, everytime a split is visited (one node having two or more outgoing edges), the current state of that path is added to unfinishedPaths. The followSubPath method consumes the unfinished paths by extending the modified DFS until all possible paths (excluding cases addressed above) are visited and completed.
        """
        while stack:
            currentNode = stack.pop()
            if type(currentNode) == list:  # handling case of second iteration on.
                currentNode = currentNode[0]

            if currentNode not in visited:  # see doc string comment for handling cases
                visited.append(currentNode)
                isLastNode = self.checkIfLastNode(currentNode=currentNode)
                if isLastNode == True:
                    yield visited
                else:
                    hasChildren, children = self.lookForChildren(
                        currentNode=currentNode
                    )
                    if hasChildren == False:
                        if self.graph[currentNode][0] not in visited:
                            stack.extend(self.graph[currentNode])
                    if hasChildren:
                        unfinishedPaths.append((visited, children))

        for path, children in unfinishedPaths:
            originalPath = list(path)
            for child in children:
                for finishedPath in self.followSubPath(
                    child, originalPath=originalPath
                ):
                    yield finishedPath

    # Input: graph (edge list)
    # Output: contiguous sequences
    def createContigs(self):
        edgesCount, startNodes = self.findStartNodes()
        incoming = sum(1 for counts in edgesCount.values() if counts[0] == 0)
        outgoing = sum(1 for counts in edgesCount.values() if counts[1] == 0)

        # file creation added for logging purposes
        edgesCountFile = os.path.join(self.logsDataDir, "edgesCount.json")
        try:
            with open(edgesCountFile, "w") as file:
                json.dump(edgesCount, file)
        except FileNotFoundError:
            print("File or directory not found")
        logging.info(f"Create Contigs: ")
        logging.info(f"\tNumber of start nodes: {incoming}")
        logging.info(f"\tNumber of end nodes: {outgoing}")

        # given the list of startNodes, explore all possible paths from each start node
        contigs = []
        for index, node in enumerate(startNodes):
            for id, path in enumerate(self.followPath(node)):
                self.allPaths.append(path)
            # logging.info(f"There were {id+1} paths for node: {index+1}")

        # for each of the paths visited, concatentate all nodes by taking the last base from the end of each node and appending it to the ongoing list of characters (initialized with the first full node)
        for path in self.allPaths:
            contig = []
            contigStr = ""
            for node in path:
                if len(contig) == 0:
                    contig.append(node)
                else:
                    contig.append(node[-1])
            contigStr = "".join(contig)
            contigs.append(contigStr)

        contigsFile = os.path.join(self.outputDataDir, "contigs.txt")
        with open(contigsFile, "w") as file:
            for contig in contigs:
                file.write(contig + "\n")

        try:
            avgLen = sum(len(contig) for contig in contigs) / len(contigs)
            logging.info(f"\tAverage contig length: {avgLen}")
            logging.info(f"\tMinimum contig length: {len(min(contigs, key=len))}")
            logging.info(f"\tMaximum contig length: {len(max(contigs, key=len))}")
            logging.info(f"\tTotal number of contigs: {[len(contigs)]}")
        except ZeroDivisionError:
            print("Length of contigs is 0, cannot calculate avg length of contig")
            print(f"len contigs: {len(contigs)}. startnodes: {startNodes}")

        return contigs
