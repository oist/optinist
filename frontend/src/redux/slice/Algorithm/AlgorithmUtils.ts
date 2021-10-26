import { OutputPath, OutputPathType } from './AlgorithmType'

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

export function isImageOutput(
  path: OutputPathType,
): path is OutputPath<'image'> {
  if (path.type === 'image') {
    return true
  } else {
    return false
  }
}

export function isPlotDataOutput(
  path: OutputPathType,
): path is OutputPath<'plotData'> {
  if (path.type === 'plotData') {
    return true
  } else {
    return false
  }
}
