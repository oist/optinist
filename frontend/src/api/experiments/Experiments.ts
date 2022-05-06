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
      position: {
        x: number
        y: number
      }
      success: string
      unique_id: string
    }
  }
  name: string
  success: string
  timestamp: string
  unique_id: string
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

export async function downloadExperimentNwbApi(uid: string) {
  const response = await axios.get(
    `${BASE_URL}/experiments/download/nwb/${uid}`,
    {
      responseType: 'blob',
    },
  )
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
