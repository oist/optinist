import axios from 'axios'

import { BASE_URL } from 'const/API'

import type { Edge, Node } from 'react-flow-renderer'
import type {
  AlgorithmNodeData,
  InputNodeData,
} from 'store/slice/FlowElement/FlowElementType'
import type { FILE_TYPE } from 'store/slice/InputNode/InputNodeType'
import type { ParamMap } from 'store/utils/param/ParamType'

export type RunPostData = {
  nodeList: Node<NodePostDataType>[]
  edgeList: Edge[]
  nwbParam: ParamMap
  snakemakeParam: ParamMap
}

export type NodePostDataType = AlgorithmNodePostData | InputNodePostData

export interface InputNodePostData extends InputNodeData {
  path: string
  fileType: FILE_TYPE
  param?: {
    [key: string]: unknown
  }
  hdf5Path?: string
}

export interface AlgorithmNodePostData extends AlgorithmNodeData {
  path: string
  param: ParamMap
}

export async function runApi(data: RunPostData): Promise<string> {
  const response = await axios.post(`${BASE_URL}/run`, data)
  return response.data
}

export async function runByUidApi(
  uid: string,
  data: RunPostData,
): Promise<string> {
  const response = await axios.post(`${BASE_URL}/run/${uid}`, data)
  return response.data
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
