export const HDF5_SLICE_NAME = 'hdf5'

export type TreeNodeTypeDTO = DirNodeDTO | FileNodeDTO

// export type ParamType = ParamParent | ParamChild

// export type ParamParent = {
//   type: 'parent'
//   children: {
//     [key: string]: ParamType
//   }
// }

// export type ParamChild = {
//   type: 'child'
//   value: unknown
//   path: string
// }

// export type ParamDTO = {
//   [key: string]: unknown
// }

export interface DirNodeDTO {
  isdir: true
  nodes: {
    [key: string]: TreeNodeTypeDTO
  }
}

export interface FileNodeDTO {
  isdir: false
}

// export type TreeNodeType = DirNode | FileNode

// export interface NodeBase {
//   name: string
//   isDir: boolean
// }

// export interface DirNode extends NodeBase {
//   isDir: true
//   nodes: TreeNodeType[]
// }

// export interface FileNode extends NodeBase {
//   isDir: false
// }

// export interface HDF5Tree {
//   isLatest: boolean
//   isLoading: boolean
//   tree: TreeNodeType[]
// }

export const HDF5_TYPE_SET = {
  ALL: 'all',
} as const

export type HDF5_TYPE = typeof HDF5_TYPE_SET[keyof typeof HDF5_TYPE_SET]
