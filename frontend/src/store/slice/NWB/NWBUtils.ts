import {
  NWBListDTO,
  NWBListType,
  NWBChild,
  NWBParent,
  NWBNodeType,
} from './NWBType'

const PATH_SEPARATOR = '/'

export function convertToNWBListType(dto: NWBListDTO, nameList?: string[]) {
  const nwbList: NWBListType = {}
  Object.entries(dto).forEach(([name, value]) => {
    const nList = nameList ?? []
    if (Object.prototype.hasOwnProperty.call(value, 'children')) {
      nList.push(name)
      nwbList[name] = {
        type: 'parent',
        children: convertToNWBListType(
          (
            value as {
              children: NWBListDTO
            }
          ).children as NWBListDTO,
          nList,
        ),
      }
    } else {
      nwbList[name] = {
        type: 'child',
        value,
        path: nList.concat([name]).join(PATH_SEPARATOR),
      }
    }
  })
  return nwbList
}

export function isNWBChild(nwbNode: NWBNodeType): nwbNode is NWBChild {
  return nwbNode.type === 'child'
}

export function isNWBParent(ndwNode: NWBNodeType): ndwNode is NWBParent {
  return ndwNode.type === 'parent'
}

export function getNWBChild(
  nwbNodeList: NWBListType,
  paramPath: string,
): NWBChild | null {
  let result: NWBChild | null = null
  for (const node of Object.values(nwbNodeList)) {
    if (isNWBChild(node)) {
      if (node.path === paramPath) {
        result = node
      }
    } else {
      result = getNWBChild(node.children, paramPath)
    }
    if (result != null) {
      break
    }
  }
  return result
}
