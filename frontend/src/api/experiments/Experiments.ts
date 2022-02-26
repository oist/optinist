import axios from 'axios'

import { BASE_URL } from 'const/API'

export type ExperimentsDTO = {
  [uid: string]: ExperimentDTO
}

export type ExperimentDTO = {
  function: {
    [name: string]: {
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

export async function getExperiments(): Promise<ExperimentsDTO> {
  const response = await axios.get(`${BASE_URL}/experiments`)
  return response.data
}

// export async function getExperimentByUid(uid: string): Promise<{
//   // todo
// }> {
//   const response = await axios.get(`${BASE_URL}/experiments/${uid}`)
//   return response.data
// }

export async function deleteExperimentByUid(uid: string): Promise<boolean> {
  const response = await axios.delete(`${BASE_URL}/experiments/${uid}`)
  return response.data
}
