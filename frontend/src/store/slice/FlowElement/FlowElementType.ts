import { Transform, Node, Edge } from 'react-flow-renderer'

export const FLOW_ELEMENT_SLICE_NAME = 'flowElement'

export const NODE_TYPE_SET = {
  INPUT: 'input',
  ALGORITHM: 'algorithm',
} as const

export type NODE_TYPE = (typeof NODE_TYPE_SET)[keyof typeof NODE_TYPE_SET]

export type NodeData = AlgorithmNodeData | InputNodeData

interface NodeDataBaseType<T extends NODE_TYPE> {
  label: string
  type: T
}

export interface InputNodeData
  extends NodeDataBaseType<typeof NODE_TYPE_SET.INPUT> {}

export interface AlgorithmNodeData
  extends NodeDataBaseType<typeof NODE_TYPE_SET.ALGORITHM> {}

export interface ElementCoord {
  x: number
  y: number
}

export interface FlowElement {
  flowNodes: Node<NodeData>[]
  flowEdges: Edge[]
  flowPosition: Transform
  elementCoord: ElementCoord
}
