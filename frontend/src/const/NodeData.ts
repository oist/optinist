interface BaseNodeData {
  label: string
  type: NODE_DATA_TYPE
}

export type NodeData = InputNodeData | OutPutNodeData | AlgoNodeData

export const NODE_DATA_TYPE_SET = {
  DATA: 'data',
  OUTPUT: 'output',
  ALGO: 'algo',
} as const

export type NODE_DATA_TYPE =
  typeof NODE_DATA_TYPE_SET[keyof typeof NODE_DATA_TYPE_SET]

export interface InputNodeData extends BaseNodeData {
  path?: string
  type: 'data'
}

export interface OutPutNodeData extends BaseNodeData {
  type: 'output'
}

export interface AlgoNodeData extends BaseNodeData {
  type: 'algo'
  param?: AlgoParam
}

export type AlgoParam = {
  [name: string]: unknown
}
