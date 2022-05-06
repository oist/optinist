import { AlgoListDTO, AlgorithmInfo } from 'api/algolist/AlgoList'
import {
  AlgorithmNodeType,
  AlgorithmChild,
  AlgorithmParent,
  AlgorithmListTree,
} from './AlgorithmListType'

export function isAlgoChild(
  algoNode: AlgorithmNodeType,
): algoNode is AlgorithmChild {
  return algoNode.type === 'child'
}

export function isAlgoParent(
  algoNode: AlgorithmNodeType,
): algoNode is AlgorithmParent {
  return algoNode.type === 'parent'
}

export function getAlgoChild(
  algoList: AlgorithmListTree,
  algoName: string,
): AlgorithmChild | null {
  let result: AlgorithmChild | null = null
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
  const algoList: AlgorithmListTree = {}
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
      const dto = value as {
        args: AlgorithmInfo[]
        returns: AlgorithmInfo[]
        path: string
      }
      algoList[name] = {
        type: 'child',
        functionPath: dto.path,
        args: dto.args,
        returns: dto.returns,
      }
    }
  })
  return algoList
}
