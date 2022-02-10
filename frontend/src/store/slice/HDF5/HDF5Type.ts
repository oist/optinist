export const HDF5_SLICE_NAME = 'hdf5'

export type TreeNodeTypeDTO = DirNodeDTO | FileNodeDTO

export interface NodeBaseDTO {
  path: string
  name: string
  isdir: boolean
}

export interface DirNodeDTO extends NodeBaseDTO {
  isdir: true
  nodes: TreeNodeTypeDTO[]
}

export interface FileNodeDTO extends NodeBaseDTO {
  isdir: false
}

export type TreeNodeType = DirNode | FileNode

export interface NodeBase {
  name: string
  isDir: boolean
}

export interface DirNode extends NodeBase {
  isDir: true
  nodes: TreeNodeType[]
}

export interface FileNode extends NodeBase {
  isDir: false
}

export interface HDF5Tree {
  isLatest: boolean
  isLoading: boolean
  tree: TreeNodeType[]
}

export const HDF5_TYPE_SET = {
  ALL: 'all',
} as const

export type HDF5_TYPE = typeof HDF5_TYPE_SET[keyof typeof HDF5_TYPE_SET]
