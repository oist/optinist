import { ParamMap } from 'store/utils/param/ParamType'

export const ALGORITHM_NODE_SLICE_NAME = 'algorithmNode'

export type AlgorithmNode = {
  [nodeId: string]: AlgrithmNodeType
}

type AlgrithmNodeType = {
  functionPath: string
  name: string
  params: ParamMap | null
  isUpdated: boolean
}
