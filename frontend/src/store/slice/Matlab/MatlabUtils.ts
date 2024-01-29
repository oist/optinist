import { MatlabTreeDTO } from "api/matlab/Matlab"
import { MatlabTreeNodeType } from "store/slice/Matlab/MatlabType"

export function convertToTreeNodeType(
  dto: MatlabTreeDTO[],
): MatlabTreeNodeType[] {
  return dto.map((node) =>
    node.isDir
      ? {
          name: node.name,
          isDir: true,
          nodes: convertToTreeNodeType(node.nodes),
          path: node.path,
          dataType: node.dataType,
        }
      : {
          name: node.name,
          isDir: false,
          shape: node.shape,
          path: node.path,
          nbytes: node.nbytes,
          dataType: node.dataType,
        },
  )
}
