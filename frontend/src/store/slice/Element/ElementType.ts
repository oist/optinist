import { Elements } from 'react-flow-renderer'
import { NodeData } from 'const/NodeData'

export const ELEMENT_SLICE_NAME = 'element'

export interface Element {
  flowElements: Elements<NodeData>
  runStatus: RUN_STATUS_TYPE
  runMessage?: string
}

export const RUN_STATUS = {
  RUNNING: 'running',
  SUCCESS: 'success',
  FAILED: 'failed',
  STOPPED: 'stopped',
} as const

export type RUN_STATUS_TYPE = typeof RUN_STATUS[keyof typeof RUN_STATUS]
