import type { Edge, Node } from "reactflow"

import { BASE_URL } from "const/API"
import type {
  AlgorithmNodeData,
  InputNodeData,
} from "store/slice/FlowElement/FlowElementType"
import type { FILE_TYPE } from "store/slice/InputNode/InputNodeType"
import axios from "utils/axios"
import type { ParamMap } from "utils/param/ParamType"

export type RunPostData = {
  name: string
  nodeDict: NodeDict
  edgeDict: EdgeDict
  nwbParam: ParamMap
  snakemakeParam: ParamMap
  forceRunList: { nodeId: string; name: string }[]
}

export type NodeDict = {
  [nodeId: string]: Node<NodePostDataType>
}

export type EdgeDict = {
  [nodeId: string]: Edge
}

export type NodePostDataType = AlgorithmNodePostData | InputNodePostData

export interface InputNodePostData extends InputNodeData {
  path: string | string[]
  fileType: FILE_TYPE
  param?: {
    [key: string]: unknown
  }
  hdf5Path?: string
  matPath?: string
}

export interface AlgorithmNodePostData extends AlgorithmNodeData {
  path: string
  param: ParamMap
}

export async function runApi(
  workspaceId: number,
  data: RunPostData,
): Promise<string> {
  const response = await axios.post(`${BASE_URL}/run/${workspaceId}`, data)
  return response.data
}

export async function runByUidApi(
  workspaceId: number,
  uid: string,
  data: Omit<RunPostData, "name">,
): Promise<string> {
  const response = await axios.post(
    `${BASE_URL}/run/${workspaceId}/${uid}`,
    data,
  )
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
  workspaceId: number
  uid: string
  pendingNodeIdList: string[]
}): Promise<RunResultDTO> {
  const { workspaceId, uid, pendingNodeIdList } = data
  const response = await axios.post(
    `${BASE_URL}/run/result/${workspaceId}/${uid}`,
    {
      pendingNodeIdList,
    },
  )
  return response.data
}

export async function cancelResultApi(data: {
  workspaceId: number
  uid: string
}): Promise<RunResultDTO> {
  const { workspaceId, uid } = data
  const response = await axios.post(
    `${BASE_URL}/run/cancel/${workspaceId}/${uid}`,
  )
  return response.data
}
