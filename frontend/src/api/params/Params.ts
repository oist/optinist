import axios from 'utils/axios'

import { BASE_URL } from 'const/API'
import { ParamDTO } from 'utils/param/ParamType'

export async function getAlgoParamsApi(algoName: string): Promise<ParamDTO> {
  const response = await axios.get(`${BASE_URL}/params/${algoName}`)
  return response.data
}
