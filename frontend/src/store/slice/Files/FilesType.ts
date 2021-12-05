export const FILES_SLICE_NAME = 'files'

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
  path: string
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

export interface Files {
  isLatest: boolean
  isLoading: boolean
  tree: TreeNodeType[]
}
