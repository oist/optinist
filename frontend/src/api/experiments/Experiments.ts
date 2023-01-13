import axios from 'axios'

import { BASE_URL } from 'const/API'
import { RunPostData } from 'api/run/Run'

export type ExperimentsDTO = {
  [uid: string]: ExperimentDTO
}

export type ExperimentDTO = {
  function: {
    [uid: string]: {
      name: string
      success: string
      unique_id: string
      hasNWB: boolean
    }
  }
  name: string
  success: string
  timestamp: string
  unique_id: string
  hasNWB: boolean
}

export async function getExperimentsApi(): Promise<ExperimentsDTO> {
  const response = await axios.get(`${BASE_URL}/experiments`)
  return response.data
}

export async function deleteExperimentByUidApi(uid: string): Promise<boolean> {
  const response = await axios.delete(`${BASE_URL}/experiments/${uid}`)
  return response.data
}

export async function deleteExperimentByListApi(
  uidList: Array<string>,
): Promise<boolean> {
  const response = await axios.post(`${BASE_URL}/experiments/delete`, {
    uidList,
  })
  return response.data
}

export async function importExperimentByUidApi(
  uid: string,
): Promise<RunPostData> {
  const response = await axios.get(`${BASE_URL}/experiments/import/${uid}`)
  return response.data
}

export async function downloadExperimentNwbApi(uid: string, nodeId?: string) {
  const path =
    nodeId != null
      ? `${BASE_URL}/experiments/download/nwb/${uid}/${nodeId}`
      : `${BASE_URL}/experiments/download/nwb/${uid}`
  const response = await axios.get(path, {
    responseType: 'blob',
  })
  return response.data
}

export async function downloadExperimentConfigApi(uid: string) {
  const response = await axios.get(
    `${BASE_URL}/experiments/download/config/${uid}`,
    {
      responseType: 'blob',
    },
  )
  return response.data
}

export async function renameExperiment(unique_id: string, new_name: string) {
  const response = await axios.patch(
    `${BASE_URL}/experiments/${unique_id}/rename`,
    {
      new_name,
    },
  )
  return response.data
}
