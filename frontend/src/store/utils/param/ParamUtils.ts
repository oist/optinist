import {
  ParamDTO,
  ParamChild,
  ParamParent,
  ParamType,
  ParamMap,
} from './ParamType'

export function getChildParam(
  path: string,
  ParamMap: ParamMap,
): ParamChild | null {
  let result: ParamChild | null = null
  for (const node of Object.values(ParamMap)) {
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

export function convertToParamMap(dto: ParamDTO, keyList?: string[]): ParamMap {
  const ParamMap: ParamMap = {}
  Object.entries(dto).forEach(([name, value]) => {
    const kList = keyList ?? []
    if (isDictObject(value)) {
      kList.push(name)
      ParamMap[name] = {
        type: 'parent',
        children: convertToParamMap(value, kList),
      }
    } else {
      ParamMap[name] = {
        type: 'child',
        value,
        path: kList.concat([name]).join(PATH_SEPARATOR),
      }
    }
  })
  return ParamMap
}
