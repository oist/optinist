import { OutputPathsDTO } from "api/run/Run"
import { BASE_URL } from "const/API"
import { EXPERIMENTS_STATUS } from "store/slice/Experiments/ExperimentsType"
import axios from "utils/axios"

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

type NWBType = {
  imaging_plane: {
    imaging_rate: number
  }
}

export type ExperimentDTO = {
  function: FunctionsDTO
  name: string
  success?: EXPERIMENTS_STATUS
  started_at: string
  finished_at?: string
  workspace_id: number
  unique_id: string
  hasNWB: boolean
  nwb: NWBType
}

export async function getExperimentsApi(
  workspaceId: number,
): Promise<ExperimentsDTO> {
  const response = await axios.get(`${BASE_URL}/experiments/${workspaceId}`)
  return response.data
}

export async function deleteExperimentByUidApi(
  workspaceId: number,
  uid: string,
): Promise<boolean> {
  const response = await axios.delete(
    `${BASE_URL}/experiments/${workspaceId}/${uid}`,
  )
  return response.data
}

export async function deleteExperimentByListApi(
  workspaceId: number,
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

export async function downloadExperimentNwbApi(
  workspaceId: number,
  uid: string,
  nodeId?: string,
) {
  const path =
    nodeId != null
      ? `${BASE_URL}/experiments/download/nwb/${workspaceId}/${uid}/${nodeId}`
      : `${BASE_URL}/experiments/download/nwb/${workspaceId}/${uid}`
  const response = await axios.get(path, {
    responseType: "blob",
  })
  return response.data
}

export async function downloadExperimentConfigApi(
  workspaceId: number,
  uid: string,
) {
  const response = await axios.get(
    `${BASE_URL}/experiments/download/config/${workspaceId}/${uid}`,
    {
      responseType: "blob",
    },
  )
  return response.data
}

export async function renameExperimentApi(
  workspaceId: number,
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
