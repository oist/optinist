import { AlgoParam } from 'const/NodeData'

export const ALGORITHM_SLICE_NAME = 'algorithm'

export type Algorithm = {
  currentAlgoId: string
  algoNodeMap: {
    [id: string]: {
      name: string
      param?: AlgoParam
      output?: OutputPaths
      selectedOutputKey?: string // key of OutputPaths
    }
  }
  algoList: AlgoListType
}

export type AlgoListType = {
  [algoName: string]: AlgoNodeType
}

export type AlgoNodeType = AlgoChild | AlgoParent
export type AlgoChild = {
  type: 'child'
  args: string[]
}
export type AlgoParent = {
  type: 'parent'
  children: {
    [name: string]: AlgoNodeType
  }
}

export type OutputPaths = {
  [key: string]: OutputPathType
}
export const OUTPUT_TYPE_SET = {
  IMAGE: 'Image',
  TIME_SERIES: 'TimeSeries',
  HEAT_MAP: 'HeatMap',
} as const

export type OUTPUT_TYPE = typeof OUTPUT_TYPE_SET[keyof typeof OUTPUT_TYPE_SET]

export type OutputPathType =
  | OutputPath<typeof OUTPUT_TYPE_SET.IMAGE>
  | OutputPath<typeof OUTPUT_TYPE_SET.TIME_SERIES>
  | OutputPath<typeof OUTPUT_TYPE_SET.HEAT_MAP>

export interface OutputPath<T extends OUTPUT_TYPE> {
  type: T
  path: T extends typeof OUTPUT_TYPE_SET.IMAGE
    ? ImagePathType
    : T extends typeof OUTPUT_TYPE_SET.TIME_SERIES
    ? TimeSeriesPathType
    : HeatMapPathType
}

type Path = { value: string }

type ImagePathType = Path & { maxIndex: number }

type TimeSeriesPathType = Path

type HeatMapPathType = Path
