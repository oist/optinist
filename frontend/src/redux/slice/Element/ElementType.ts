import { Elements } from 'react-flow-renderer'

export const ELEMENT_SLICE_NAME = 'element'

export type Param = {
  [name: string]: unknown
}

export type Algorithm = {
  [id: string]: {
    name: string
    param: Param
  }
}

export interface Element {
  flowElements: Elements<NodeDataType>
  clickedNodeId: string | null
  currentAlgoId: string
  algoParams: Algorithm
}

interface NodeData {
  label: string
  type: NodeType
}

export type NodeDataType = InputNodeData | OutPutNodeData | AlgoNodeData

export type NodeType = 'input' | 'output' | 'algo'

export interface InputNodeData extends NodeData {
  path?: string
  type: 'input'
}

export interface OutPutNodeData extends NodeData {
  type: 'output'
}

export interface AlgoNodeData extends NodeData {
  type: 'algo'
}
