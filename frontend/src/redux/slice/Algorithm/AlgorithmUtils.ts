import { OutputPath, OutputPathType, OUTPUT_TYPE_SET } from './AlgorithmType'

export function isImageOutput(
  path: OutputPathType,
): path is OutputPath<typeof OUTPUT_TYPE_SET.IMAGE> {
  if (path.type === OUTPUT_TYPE_SET.IMAGE) {
    return true
  } else {
    return false
  }
}

export function isTimeSeriesOutput(
  path: OutputPathType,
): path is OutputPath<typeof OUTPUT_TYPE_SET.TIME_SERIES> {
  if (path.type === OUTPUT_TYPE_SET.TIME_SERIES) {
    return true
  } else {
    return false
  }
}

export function isHeatMapOutput(
  path: OutputPathType,
): path is OutputPath<typeof OUTPUT_TYPE_SET.HEAT_MAP> {
  if (path.type === OUTPUT_TYPE_SET.HEAT_MAP) {
    return true
  } else {
    return false
  }
}
