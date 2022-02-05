import { Elements, FlowTransform } from 'react-flow-renderer'

export const FLOW_ELEMENT_SLICE_NAME = 'flowElement'

export const NODE_TYPE_SET = {
  INPUT: 'input',
  ALGORITHM: 'algorithm',
} as const

export type NODE_TYPE = typeof NODE_TYPE_SET[keyof typeof NODE_TYPE_SET]

export interface NodeData {
  label: string
  type: NODE_TYPE
}

export interface FlowElement {
  flowElements: Elements<NodeData>
  flowPosition: FlowTransform
}
