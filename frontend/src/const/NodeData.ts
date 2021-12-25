interface BaseNodeData {
  label: string
  type: NODE_DATA_TYPE
}

export type NodeData = ImageNodeData | AlgoNodeData | CsvNodeData

export const NODE_DATA_TYPE_SET = {
  IMAGE: 'image',
  CSV: 'csv',
  ALGO: 'algo',
} as const

export type NODE_DATA_TYPE =
  typeof NODE_DATA_TYPE_SET[keyof typeof NODE_DATA_TYPE_SET]

export interface ImageNodeData extends BaseNodeData {
  path?: string
  type: typeof NODE_DATA_TYPE_SET.IMAGE
}

export interface CsvNodeData extends BaseNodeData {
  path?: string
  type: typeof NODE_DATA_TYPE_SET.CSV
}

export interface AlgoNodeData extends BaseNodeData {
  type: typeof NODE_DATA_TYPE_SET.ALGO
  param?: AlgoParam
  path: string
}

export type AlgoParam = {
  [name: string]: unknown
}
