export const ALGORITHM_NODE_SLICE_NAME = 'algorithmNode'

export type AlgorithmNode = {
  [nodeId: string]: AlgrithmNodeType
}

export type AlgorithmParam = {
  [key: string]: unknown
}

type AlgrithmNodeType = {
  functionPath: string
  name: string
  params: AlgorithmParam | null
  selectedOutputKey: string | null
}
