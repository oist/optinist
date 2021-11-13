import {
  AlgoChild,
  AlgoListType,
  AlgoNodeType,
  AlgoParent,
  OutputPath,
  OutputPathType,
  OUTPUT_TYPE_SET,
} from './AlgorithmType'

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

export function isAlgoChild(algoNode: AlgoNodeType): algoNode is AlgoChild {
  return algoNode.type === 'child'
}

export function isAlgoParent(algoNode: AlgoNodeType): algoNode is AlgoParent {
  return algoNode.type === 'parent'
}

export function getAlgoChild(
  algoList: AlgoListType,
  algoName: string,
): AlgoChild | null {
  let result: AlgoChild | null = null
  Object.entries(algoList).forEach(([name, node]) => {
    if (isAlgoChild(node)) {
      if (name === algoName) {
        result = node
      }
    } else {
      result = getAlgoChild(node.children, name)
    }
  })
  return result
}
