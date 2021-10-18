import { AlgoParam } from 'const/NodeData'
import { AlgoOutputDataDTO } from './AlgorithmUtils'

export const ALGORITHM_SLICE_NAME = 'algorithm'

export type Algorithm = {
  currentAlgoId: string
  algoMap: {
    [id: string]: {
      name: string
      param?: AlgoParam
      output?: OutputPaths
      selectedPath: {
        value: string | null
        isImage: boolean
      } | null
    }
  }
  plotDataMap: {
    [id: string]: AlgoOutputDataDTO // todo 後で型を検討
  }
}

export type OutputData = {
  xLabels: string[]
  yValues: number[]
  legends: string[]
}

export type OutputPaths = {
  images?: {
    maxIndex: number
    path: string
  }
  fluo?: string
}
