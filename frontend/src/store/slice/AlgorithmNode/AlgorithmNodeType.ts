import { ParamMap } from "utils/param/ParamType"

export const ALGORITHM_NODE_SLICE_NAME = "algorithmNode"

export type AlgorithmNode = {
  [nodeId: string]: AlgorithmNodeType
}

type AlgorithmNodeType = {
  functionPath: string
  name: string
  params: ParamMap | null
  originalValue: unknown
  isUpdate: boolean
}
