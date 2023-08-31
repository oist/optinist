import { Edge, Node, Viewport } from 'reactflow'

export const FLOW_ELEMENT_SLICE_NAME = 'flowElement'

export const NODE_TYPE_SET = {
  INPUT: 'input',
  ALGORITHM: 'algorithm',
} as const

export type NODE_TYPE = typeof NODE_TYPE_SET[keyof typeof NODE_TYPE_SET]

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
  flowElements: (Node<NodeData> | Edge<NodeData>)[]
  flowPosition: Viewport
  elementCoord: ElementCoord
}
