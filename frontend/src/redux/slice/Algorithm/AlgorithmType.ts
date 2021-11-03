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
    // idの方はOutputDataIdを使いたいけど、テンプレートリテラルはオブジェクトのkeyとして使えない...
    [id: string]: AlgoOutputDataDTO // todo 後で型を検討
  }
}

/**
 * ${nodeId}/${outputKey}
 */
export type OutputDataId = `${string}/${string}`

export type OutputData = {
  xLabels: string[]
  yValues: number[]
  legends: string[]
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

// type OutputType = 'image' | 'timeseries' | 'heatMap'

type Path = { value: string }

type ImagePathType = Path & { maxIndex: number }

type TimeSeriesPathType = Path

type HeatMapPathType = Path
