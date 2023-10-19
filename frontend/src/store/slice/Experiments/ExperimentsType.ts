export const EXPERIMENTS_SLICE_NAME = "experiments"

export type Experiments =
  | {
      status: "fulfilled"
      experimentList: ExperimentListType
    }
  | {
      status: "uninitialized"
    }
  | {
      status: "pending"
    }
  | {
      status: "error"
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
  status?: EXPERIMENTS_STATUS
  name: string
  startedAt: string
  finishedAt?: string
  hasNWB: boolean
  frameRate?: number
}

export type ExperimentFunction = {
  name: string
  nodeId: string
  status: EXPERIMENTS_STATUS
  hasNWB: boolean
  message?: string
}

export type EXPERIMENTS_STATUS = "success" | "error" | "running"

export interface ExperimentSortKeys {
  uid: string
  name: string
  startedAt: string
}
