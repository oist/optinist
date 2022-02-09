import { TreeNodeTypeDTO, TreeNodeType } from './hdf5TreeType'

export function convertToTreeNodeType(dto: TreeNodeTypeDTO[]): TreeNodeType[] {
  return dto.map((node) =>
    node.isdir
      ? {
          path: node.path,
          name: node.name,
          isDir: true,
          nodes: convertToTreeNodeType(node.nodes),
        }
      : {
          path: node.path,
          name: node.name,
          isDir: false,
        },
  )
}
