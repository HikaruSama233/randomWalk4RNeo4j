library(RNeo4j)

importGraphData = function(fileName, graph){
  clear(graph, input = FALSE)
  dt = read.table(fileName, header = FALSE, sep = "\t")
  names(dt) = c("FromNodeId", "ToNodeId")
  csvFile = paste(fileName, ".csv", sep = "")
  write.csv(dt, , csvFile, row.names = FALSE)
  addConstraint(graph, "Node", "id")
  # use cypher to load nodes and links from CSV file
  query = sprintf("
  USING PERIODIC COMMIT 50000
  LOAD CSV WITH HEADERS FROM 'file:///%s' AS df
  MERGE (from:Node {id: df.FromNodeId})
  MERGE (to:Node {id: df.ToNodeId})
  CREATE (from)-[fof:FRIENDORFOE]->(to)
  ", csvFile)
  cypher(graph, query) 
  
  summary(graph)
  print(getIndex(graph))
  print(getConstraint(graph))
  graph
}

isVisited = function(nextV) {
  # check if a node is visited
  flag = FALSE
  if ("visited" %in% names(nextV)) {
    if (as.integer(nextV$visited) > 0) {
      flag = TRUE
    }
  }
  else {
    updateProp(nextV, visited = 0)
  }
  flag
}

removeVisitedRecord = function(graph){
  nodeList = getAllVisitedV(graph)
  rm = lapply(nodeList, function(x) updateProp(x, visited = 0))
}

getAllVisitedV = function(graph){
  query = "
  MATCH (n:Node)
  WHERE n.visited = 1
  RETURN n
  "
  nodeList = getNodes(graph, query)
  nodeList
}
updateSRWpasses = function(nextV){
  if (isVisited(nextV)){
    newSRWpasses = as.integer(nextV$SRWpasses) + 1
    updateProp(nextV, SRWpasses = newSRWpasses)
  }
  else {
    updateProp(nextV, SRWpasses = 1, visited = 1)
  }
}

singleRandomWalk = function(graph, startVid, steps){
  # startVid: the id of start vertex
  # steps: the number of steps to make
  startVid = as.character(startVid)
  startV = getSingleNode(graph,
                     "
                     MATCH (n:Node {id:{svid}})
                     RETURN n",
                    svid = startVid)
  updateSRWpasses(startV)
  outPath = rep(NA, (steps + 1))
  outPath[1] = startVid
  for (i in 1:steps){
    query = "
      MATCH (from: Node {id: {fid}})-[fof:FRIENDORFOE]->(to: Node)
      RETURN to
      "
    nodeList = getNodes(graph, query, fid = startVid)
    
    # randomly select next node
    nextV = nodeList[[sample(c(1:length(nodeList)), size = 1)]]
    if ((length(nodeList) == 1) & (startVid == nextV$id)){
      break
    }
    updateSRWpasses(nextV)
    startVid = nextV$id
    outPath[i + 1] = startVid
  }
  if (i < steps) {
    outPath[c(1:i)]
  }
  else{
    outPath
  }
}

multiRandomWalk = function(graph, numStart, steps){
  # numStart: The number of different start points
  
  totalNumOfNodes = cypher(graph, "MATCH (n:Node) RETURN COUNT(n)")[1, 1]
  allIds = c(0:(totalNumOfNodes - 1)) # can be changed into real ids
  multiOutPath = lapply(allIds[sample(c(1:totalNumOfNodes), numStart)],
                        function(x) singleRandomWalk(graph, x, steps))
  # multiOutPath: list of pathes
  outPassVectors = as.data.frame(matrix(0, nrow = totalNumOfNodes, ncol = numStart))
  row.names(outPassVectors) = allIds
  for (i in 1:numStart) {
    temp = table(multiOutPath[[i]])
    for (j in 1:length(temp)){
      x = names(temp)[j]
      outPassVectors[x, i] = temp[x][[1]]
    }
  }
  
  list(paths = multiOutPath, passVectors = outPassVectors)
}

inRopeLength = function(graph, nodeList, ropeLength){
  maxL = 0
  for (i in 1:(length(nodeList) - 1)){
    for (j in (i + 1):length(nodeList)){
      query = sprintf("
        MATCH (from:Node {id:'%s'}), (to:Node {id:'%s'}), 
        p = shortestPath((from)-[*]-(to))
        RETURN p
        ", nodeList[i], nodeList[j])
      p = cypherToList(graph, query)
      pathL = ifelse(length(p) == 0, 0, p[[1]]$p$length) # length of shortest path
      
      if (pathL > maxL){
        maxL = pathL
      }
    }
  }
  maxL <= ropeLength
}

multiAgentRandomWalk = function(graph, startVid, numAgents, ropeLength, steps){
  # startVid: the id of start vertex
  # numAgents: number of agents
  # ropeLength: rope length
  # steps: the number of steps to make
  startVid = as.character(startVid)
  startV = getSingleNode(graph,
                         "
                         MATCH (n:Node {id:{svid}})
                         RETURN n",
                         svid = startVid)
  pathsDf = as.data.frame(matrix(NA, ncol = numAgents, nrow = (steps + 1)))
  # pathsDF: stores paths for all agents
  pathsDf[1,] = rep(startVid, numAgents)
  currentStep = 0
  query = "
    MATCH (from: Node {id: {fid}})-[fof:FRIENDORFOE]->(to: Node)
    RETURN to
  "
  while (currentStep <= steps) {  # get paths for each agent
    if (currentStep == 0){
      neighborNL = getNodes(graph, query, fid = startVid) # neighbors of a node
      neighborNLId = sapply(neighborNL, function(x) x$id)
      repeat {
        nextNodes = sample(neighborNLId, numAgents, replace = TRUE)
        if (inRopeLength(graph, nextNodes, ropeLength)){
          currentStep = currentStep + 1
          pathsDf[(currentStep + 1),] = nextNodes
          break
        }
      }
    } else {
      neighborNL = list()
      for (agent_i in 1:numAgents) {
        Vid = pathsDf[(currentStep + 1), agent_i]
        temp = getNodes(graph, query, fid = Vid)
        neighborNL[[agent_i]] = sapply(temp, function(x) x$id)
      }
      repeat {
        nextNodes = sapply(neighborNL, function(x) sample(x, 1))
        if (inRopeLength(graph, nextNodes, ropeLength)){
          currentStep = currentStep + 1
          pathsDf[(currentStep + 1),] = nextNodes
          break
        }
      }
    }
  }
  
  totalNumOfNodes = cypher(graph, "MATCH (n:Node) RETURN COUNT(n)")[1, 1]
  outPassVectors = as.data.frame(matrix(0, nrow = totalNumOfNodes, ncol = numAgents))
  allIds = c(0:(totalNumOfNodes - 1)) # can be changed into real ids
  row.names(outPassVectors) = allIds
  outPassChain = rep(0, totalNumOfNodes)
  names(outPassChain) = allIds
  for (i in 1:numAgents) {
    temp = table(pathsDf[[i]])
    for (j in 1:length(temp)){
      x = names(temp)[j]
      outPassVectors[x, i] = temp[x][[1]]
    }
  }
  for (i in 1:nrow(pathsDf)) {
    uniqV = unique(t(pathsDf[i,])) # get visited point for chain in each step
    outPassChain[uniqV] = outPassChain[uniqV] + 1
  }
  
  list(passesOfWalkers = outPassVectors, passesOfChain = outPassChain)
}

list2graph = function(graph, nodeList){
  l = length(nodeList)
  fromNodeId = rep(NA, l*l) 
  toNodeId = rep(NA, l*l)
  temp = 0
  # in the worest case, all nodes are linked with each other
  for (i in 1:l) {
    for (j in 1:l){
      query = sprintf("
        MATCH (from:Node {id:'%s'})-[r]->(to:Node {id:'%s'})
        RETURN r
        ", nodeList[i], nodeList[j])
      p = cypherToList(graph, query)
      if (length(p) > 0) {
        temp = temp + 1
        fromNodeId[temp] = nodeList[i]
        toNodeId[temp] = nodeList[j]
      }
    }
  }
  subGraphDf = as.data.frame(cbind(fromNodeId[1:temp], toNodeId[1:temp]))
  names(subGraphDf) = c("FromNodeId", "ToNodeId")
  subGraphDf
}

#connect to graph
graph = startGraph("http://localhost:7474/db/data/",
                   username,
                   password)

clear(graph, input = FALSE)
