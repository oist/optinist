import axios from 'axios'

import { BASE_URL } from 'const/API'
import { EdgeDict, NodeDict, OutputPathsDTO, RunPostData } from 'api/run/Run'
import { EXPERIMENTS_STATUS } from 'store/slice/Experiments/ExperimentsType'

export type ExperimentsDTO = {
  [uid: string]: ExperimentDTO
}

export type FunctionsDTO = {
  [nodeId: string]: {
    name: string
    success: string
    unique_id: string
    hasNWB: boolean
    message?: string
    started_at?: string
    finished_at?: string
    outputPaths?: OutputPathsDTO
  }
}

export type ExperimentDTO = {
  function: FunctionsDTO
  name: string
  success?: EXPERIMENTS_STATUS
  started_at: string
  finished_at?: string
  workspace_id: string
  unique_id: string
  hasNWB: boolean
  edgeDict: EdgeDict
  nodeDict: NodeDict
}

export async function getExperimentsApi(
  workspaceId: string,
): Promise<ExperimentsDTO> {
  const response = await axios.get(`${BASE_URL}/experiments/${workspaceId}`)
  return response.data
}

export async function deleteExperimentByUidApi(
  workspaceId: string,
  uid: string,
): Promise<boolean> {
  const response = await axios.delete(
    `${BASE_URL}/experiments/${workspaceId}/${uid}`,
  )
  return response.data
}

export async function deleteExperimentByListApi(
  workspaceId: string,
  uidList: Array<string>,
): Promise<boolean> {
  const response = await axios.post(
    `${BASE_URL}/experiments/delete/${workspaceId}`,
    {
      uidList,
    },
  )
  return response.data
}

export async function importExperimentByUidApi(
  workspaceId: string,
  uid: string,
): Promise<RunPostData> {
  const response = await axios.get(
    `${BASE_URL}/experiments/import/${workspaceId}/${uid}`,
  )
  return response.data
}

export async function downloadExperimentNwbApi(
  workspaceId: string,
  uid: string,
  nodeId?: string,
) {
  const path =
    nodeId != null
      ? `${BASE_URL}/experiments/download/nwb/${workspaceId}/${uid}/${nodeId}`
      : `${BASE_URL}/experiments/download/nwb/${workspaceId}/${uid}`
  const response = await axios.get(path, {
    responseType: 'blob',
  })
  return response.data
}

export async function downloadExperimentConfigApi(
  workspaceId: string,
  uid: string,
) {
  const response = await axios.get(
    `${BASE_URL}/experiments/download/config/${workspaceId}/${uid}`,
    {
      responseType: 'blob',
    },
  )
  return response.data
}

export async function renameExperiment(
  workspaceId: string,
  uid: string,
  new_name: string,
) {
  const response = await axios.patch(
    `${BASE_URL}/experiments/${workspaceId}/${uid}/rename`,
    {
      new_name,
    },
  )
  return response.data
}
