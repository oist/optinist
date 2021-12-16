import {
  AlgoChild,
  AlgoInfo,
  AlgoListDTO,
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
  for (const [name, node] of Object.entries(algoList)) {
    if (isAlgoChild(node)) {
      if (name === algoName) {
        result = node
      }
    } else {
      result = getAlgoChild(node.children, algoName)
    }
    if (result != null) {
      break
    }
  }
  return result
}

export function convertToAlgoListType(dto: AlgoListDTO) {
  const algoList: AlgoListType = {}
  Object.entries(dto).forEach(([name, value]) => {
    if (Object.prototype.hasOwnProperty.call(value, 'children')) {
      algoList[name] = {
        type: 'parent',
        children: convertToAlgoListType(
          (
            value as {
              children: AlgoListDTO
            }
          ).children as AlgoListDTO,
        ),
      }
    } else {
      algoList[name] = {
        type: 'child',
        ...(value as {
          args: AlgoInfo[]
          returns: AlgoInfo[]
          path: string
        }),
      }
    }
  })
  return algoList
}
