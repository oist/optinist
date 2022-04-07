def get_network(runItem):
    nodeList = runItem.nodeList
    edgeList = runItem.edgeList

    # nodeを初期化
    nodeDict = {}
    for node in nodeList:
        nodeDict[node['id']] = node

    returnCntDict = {key: 0 for key in nodeDict.keys()}
    for edge in edgeList:
        returnCntDict[edge["source"]] += 1

    endNodeList = []
    for key, value in returnCntDict.items():
        if value == 0:
            endNodeList.append(key)

    return nodeDict, edgeList, endNodeList
