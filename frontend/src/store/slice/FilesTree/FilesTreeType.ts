export const FILES_TREE_SLICE_NAME = 'filesTree'

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

export interface FilesTree {
  [fileType: string]: {
    isLatest: boolean
    isLoading: boolean
    tree: TreeNodeType[]
  }
}

export const FILE_TYPE_SET = {
  IMAGE: 'image',
  CSV: 'csv',
  ALL: 'all',
} as const

export type FILE_TYPE = typeof FILE_TYPE_SET[keyof typeof FILE_TYPE_SET]
