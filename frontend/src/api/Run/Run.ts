import axios from 'axios'

import { BASE_URL } from 'const/API'

import { Edge, Node } from 'react-flow-renderer'
import { ParamMap } from 'store/utils/param/ParamType'

export type RunPostData = {
  nodeList: Node[]
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

export const RUN_RESULT_STATUS = {
  SUCCESS: 'success',
  ERROR: 'error',
} as const

export type RUN_RESULT_STATUS_TYPE =
  typeof RUN_RESULT_STATUS[keyof typeof RUN_RESULT_STATUS]

export type RunResultDTO = {
  [path: string]: {
    path: string
    stutus: RUN_RESULT_STATUS_TYPE
  }
}

export async function runResult(data: { uid: string }): Promise<RunResultDTO> {
  const response = await axios.post(`${BASE_URL}/run/result/${data.uid}`)
  return response.data
}
