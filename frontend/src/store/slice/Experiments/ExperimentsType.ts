export const EXPERIMENTS_SLICE_NAME = 'experiments'

export type Experiments =
  | {
      status: 'fulfilled'
      experimentList: ExperimentListType
    }
  | {
      status: 'uninitialized'
    }
  | {
      status: 'pending'
    }
  | {
      status: 'error'
      message?: string
    }

export type ExperimentListType = {
  [uid: string]: ExperimentType
}

export type ExperimentType = {
  uid: string
  functions: {
    [nodeId: string]: ExperimentFunction
  }
  name: string
  timestamp: string
}

export type ExperimentFunction = {
  name: string
  nodeId: string
  postion: {
    x: number
    y: number
  }
  status: EXPERIMENTS_STATUS
}

export type EXPERIMENTS_STATUS = 'success' | 'error' | 'running'

export interface Experiment {
  uid: string
  name: string
  timestamp: string
}
