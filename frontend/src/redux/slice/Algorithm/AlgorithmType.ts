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
      selectedOutputKey?: string // key of OutputPaths
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
  [key: string]: OutputPathType
}

export type OutputPathType = OutputPath<'image'> | OutputPath<'plotData'>

export interface OutputPath<T extends OutputType> {
  type: T
  path: T extends 'image' ? ImagePathType : PlotDataPathType
}

type OutputType = 'image' | 'plotData'

type ImagePathType = { value: string; maxIndex: number }

type PlotDataPathType = { value: string }
