import { HDF5TreeDTO } from './HDF5Type'

export function convertToTreeNodeType(dto: HDF5TreeDTO[]): HDF5TreeDTO[] {
  return dto.map((node) =>
    node.isDir
      ? {
          name: node.name,
          isDir: true,
          nodes: convertToTreeNodeType(node.nodes),
          path: node.path,
        }
      : {
          name: node.name,
          isDir: false,
          shape: node.shape,
          path: node.path,
          nbytes: node.nbytes,
        },
  )
}
