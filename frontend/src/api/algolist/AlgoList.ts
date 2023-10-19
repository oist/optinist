import { BASE_URL } from "const/API"
import axios from "utils/axios"


export type AlgoListDTO = {
  [name: string]:
    | {
        args: AlgorithmInfo[]
        returns: AlgorithmInfo[]
        path: string
      }
    | { children: AlgoListDTO }
}

export type AlgorithmInfo = {
  name: string
  type: string
  isNone?: boolean
}

export async function getAlgoListApi(): Promise<AlgoListDTO> {
  const response = await axios.get(`${BASE_URL}/algolist`)
  return response.data
}
