import { RunPostData } from 'api/run/Run'
import axios from 'axios'

import { BASE_URL } from 'const/API'

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

// export async function getExperimentByUid(uid: string): Promise<{
//   // todo
// }> {
//   const response = await axios.get(`${BASE_URL}/experiments/${uid}`)
//   return response.data
// }

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

export async function downloadExperimentByUidApi(uid: string): Promise<Blob> {
  const response = await axios.get(`${BASE_URL}/experiments/download/${uid}`)
  return response.data
}
