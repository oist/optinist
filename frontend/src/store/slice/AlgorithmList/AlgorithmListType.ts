import { AlgorithmInfo } from 'api/algolist/AlgoList'

export const ALGORITHM_LIST_SLICE_NAME = 'algorithmList'

export type AlgorithmListType = {
  isLatest: boolean
  tree: AlgorithmListTree
}

export type AlgorithmListTree = {
  [algoName: string]: AlgorithmNodeType
}

export type AlgorithmNodeType = AlgorithmChild | AlgorithmParent
export type AlgorithmChild = {
  type: 'child'
  args: AlgorithmInfo[]
  returns: AlgorithmInfo[]
  functionPath: string
}
export type AlgorithmParent = {
  type: 'parent'
  children: {
    [name: string]: AlgorithmNodeType
  }
}
