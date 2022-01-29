export const ALGORITHM_NODE_SLICE_NAME = 'algorithmNode'

export type AlgorithmNode = {
  [nodeId: string]: AlgrithmNodeType
}

type AlgrithmNodeType = {
  functionPath: string
  name: string
  params: AlgorithmParam | null
  selectedOutputKey: string | null
}

export type AlgorithmParam = {
  [paramKey: string]: ParamType
}

export type ParamType = ParamParent | ParamChild

export type ParamParent = {
  type: 'parent'
  children: {
    [key: string]: ParamType
  }
}

export type ParamChild = {
  type: 'child'
  value: unknown
  path: string
}

export type ParamDTO = {
  [key: string]: unknown
}
