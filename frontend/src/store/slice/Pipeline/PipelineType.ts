import { RunPostData } from 'api/run/Run'
import { DATA_TYPE } from '../DisplayData/DisplayDataType'

export const PIPELINE_SLICE_NAME = 'pipeline'

export type Pipeline = {
  run: PipelineType
  currentPipeline?: {
    uid: string
  }
}

export const RUN_STATUS = {
  /**
   * 未実行
   */
  START_UNINITIALIZED: 'StartUninitialized',
  /**
   * 開始リクエスト中
   */
  START_PENDING: 'StartPending',
  /**
   * 開始できなかった
   */
  START_ERROR: 'StartError',
  /**
   * 正常に開始した
   */
  START_SUCCESS: 'StartSuccess',
  /**
   * 全て終了
   */
  FINISHED: 'Finished',
  /**
   * 途中終了
   */
  ABORTED: 'Aborted',
  /**
   * キャンセルされた
   */
  CANCELED: 'Canceled',
} as const

export type RUN_STATUS_TYPE = typeof RUN_STATUS[keyof typeof RUN_STATUS]

export type PipelineType =
  | StartedPipeline
  | {
      status:
        | typeof RUN_STATUS.START_ERROR
        | typeof RUN_STATUS.START_PENDING
        | typeof RUN_STATUS.START_UNINITIALIZED
    }

export type StartedPipeline = {
  uid: string
  status:
    | typeof RUN_STATUS.START_SUCCESS
    | typeof RUN_STATUS.FINISHED
    | typeof RUN_STATUS.ABORTED
    | typeof RUN_STATUS.CANCELED
  runPostData: RunPostData
  runResult: RunResult
}

export type RunResult = {
  [nodeId: string]: NodeResult
}

export type NodeResult = NodeResultPending | NodeResultSuccess | NodeResultError

export const NODE_RESULT_STATUS = {
  SUCCESS: 'success',
  ERROR: 'error',
  PENDING: 'pending',
} as const

export type NODE_RESULT_STATUS_TYPE =
  typeof NODE_RESULT_STATUS[keyof typeof NODE_RESULT_STATUS]

export interface NodeResultBaseType<T extends NODE_RESULT_STATUS_TYPE> {
  status: T
  message?: string
  name: string
}

export interface NodeResultPending
  extends NodeResultBaseType<typeof NODE_RESULT_STATUS.PENDING> {}

export interface NodeResultSuccess
  extends NodeResultBaseType<typeof NODE_RESULT_STATUS.SUCCESS> {
  outputPaths: OutputPaths
}

export type OutputPaths = {
  [outputKey: string]: {
    path: string
    type: DATA_TYPE
  }
}

export interface NodeResultError
  extends NodeResultBaseType<typeof NODE_RESULT_STATUS.ERROR> {
  message: string
}
