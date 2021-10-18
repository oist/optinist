import { AlgoParam } from 'const/NodeData'

export const ALGORITHM_SLICE_NAME = 'algorithm'

export type Algorithm = {
  currentAlgoId: string
  algoMap: {
    [id: string]: {
      name: string
      param: AlgoParam
      output?: OutputPaths
    }
  }
}

export type AlgoOutput = {
  data: OutputData[]
}

export type OutputData = {
  x: number | string
  y: number
}

export type OutputPaths = {
  imageDir: string
}
