import { OutputData } from './AlgorithmType'

export type AlgoOutputDataDTO = {
  data: {
    [key: string]: {
      [key: number]: number
    }
  }
}

export function convertToOutputData(dto: AlgoOutputDataDTO) {
  // todo 後で変換処理を書く
  return dto
}
