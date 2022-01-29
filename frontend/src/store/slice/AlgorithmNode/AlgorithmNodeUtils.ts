import {
  ParamDTO,
  ParamChild,
  ParamParent,
  ParamType,
  AlgorithmParam,
} from './AlgorithmNodeType'

export function getChildParam(
  path: string,
  algorithmParam: AlgorithmParam,
): ParamChild | null {
  let result: ParamChild | null = null
  for (const node of Object.values(algorithmParam)) {
    if (isParamChild(node)) {
      if (node.path === path) {
        result = node
      }
    } else {
      result = getChildParam(path, node.children)
    }
    if (result != null) {
      break
    }
  }
  return result
}

export function isParamChild(param: ParamType): param is ParamChild {
  return param.type === 'child'
}

export function isParamParent(param: ParamType): param is ParamParent {
  return param.type === 'parent'
}

function isDictObject(value: unknown): value is { [key: string]: unknown } {
  return value !== null && typeof value === 'object' && !Array.isArray(value)
}

const PATH_SEPARATOR = '/'

export function convertToAlgorithmParam(
  dto: ParamDTO,
  keyList?: string[],
): AlgorithmParam {
  const algorithmParam: AlgorithmParam = {}
  Object.entries(dto).forEach(([name, value]) => {
    const kList = keyList ?? []
    if (isDictObject(value)) {
      kList.push(name)
      algorithmParam[name] = {
        type: 'parent',
        children: convertToAlgorithmParam(value, kList),
      }
    } else {
      algorithmParam[name] = {
        type: 'child',
        value,
        path: kList.concat([name]).join(PATH_SEPARATOR),
      }
    }
  })
  return algorithmParam
}
