import {
  SnakemakeListDTO,
  SnakemakeListType,
  SnakemakeChild,
  SnakemakeParent,
  SnakemakeNodeType,
} from './SnakemakeType'

const PATH_SEPARATOR = '/'

export function convertToSnakemakeListType(
  dto: SnakemakeListDTO,
  nameList?: string[],
) {
  const SnakemakeList: SnakemakeListType = {}
  Object.entries(dto).forEach(([name, value]) => {
    const nList = nameList ?? []
    if (Object.prototype.hasOwnProperty.call(value, 'children')) {
      nList.push(name)
      SnakemakeList[name] = {
        type: 'parent',
        children: convertToSnakemakeListType(
          (
            value as {
              children: SnakemakeListDTO
            }
          ).children as SnakemakeListDTO,
          nList,
        ),
      }
    } else {
      SnakemakeList[name] = {
        type: 'child',
        value,
        path: nList.concat([name]).join(PATH_SEPARATOR),
      }
    }
  })
  return SnakemakeList
}

export function isSnakemakeChild(
  SnakemakeNode: SnakemakeNodeType,
): SnakemakeNode is SnakemakeChild {
  return SnakemakeNode.type === 'child'
}

export function isSnakemakeParent(
  ndwNode: SnakemakeNodeType,
): ndwNode is SnakemakeParent {
  return ndwNode.type === 'parent'
}

export function getSnakemakeChild(
  SnakemakeNodeList: SnakemakeListType,
  paramPath: string,
): SnakemakeChild | null {
  let result: SnakemakeChild | null = null
  for (const node of Object.values(SnakemakeNodeList)) {
    if (isSnakemakeChild(node)) {
      if (node.path === paramPath) {
        result = node
      }
    } else {
      result = getSnakemakeChild(node.children, paramPath)
    }
    if (result != null) {
      break
    }
  }
  return result
}
