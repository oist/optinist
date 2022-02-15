import axios from 'axios'

import { BASE_URL } from 'const/API'

import { Edge, Node } from 'react-flow-renderer'
import { NodeData } from 'store/slice/FlowElement/FlowElementType'
import { ParamMap } from 'store/utils/param/ParamType'

export type RunPostData = {
  nodeList: Node<NodeData>[]
  edgeList: Edge[]
  nwbParam: ParamMap
  snakemakeParam: ParamMap
}

export async function run(data: {
  runData: RunPostData
  uid?: string
}): Promise<string> {
  const response = await axios.post(`${BASE_URL}/run`, data.runData)
  return response.data // uid
}

export type RunResultDTO = {
  [nodeId: string]: {
    status: string
    message: string
    name: string
    outputPaths?: OutputPathsDTO
  }
}

export type OutputPathsDTO = {
  [outputKey: string]: {
    path: string
    type: string
  }
}

export async function runResult(data: {
  uid: string
  pendingNodeIdList: string[]
}): Promise<RunResultDTO> {
  const { uid, pendingNodeIdList } = data
  const response = await axios.post(`${BASE_URL}/run/result/${uid}`, {
    pendingNodeIdList,
  })
  return response.data
}
