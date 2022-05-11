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

export function equalsParamMap(a: ParamMap, b: ParamMap) {
  if (a === b) {
    return true
  }
  const aArray = Object.keys(a)
  const bArray = Object.keys(b)
  return (
    aArray.length === bArray.length &&
    aArray.every((aKey) => {
      const aValue = a[aKey]
      const bValue = b[aKey]
      return equalsParam(aValue, bValue)
    })
  )
}

function equalsParam(a: ParamType, b: ParamType): boolean {
  if (a === b) {
    return true
  }
  if (isParamChild(a) && isParamChild(b)) {
    return equalsParamChild(a, b)
  } else if (isParamParent(a) && isParamParent(b)) {
    const aArray = Object.keys(a)
    const bArray = Object.keys(b)
    return (
      aArray.length === bArray.length &&
      aArray.every((aKey) => {
        const aValue = a.children[aKey]
        const bValue = b.children[aKey]
        return equalsParam(aValue, bValue)
      })
    )
  } else {
    return false
  }
}

function equalsParamChild(a: ParamChild, b: ParamChild) {
  return a.path === b.path && a.value === b.value
}
