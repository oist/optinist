export const FILES_TREE_SLICE_NAME = 'filesTree'

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
