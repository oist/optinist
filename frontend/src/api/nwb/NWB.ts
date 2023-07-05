import axios from 'utils/axios'

import { BASE_URL } from 'const/API'
import { ParamDTO } from 'utils/param/ParamType'

export async function getNWBParamsApi(): Promise<ParamDTO> {
  const response = await axios.get(`${BASE_URL}/nwb`)
  return response.data
}
