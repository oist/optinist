import { Elements } from 'react-flow-renderer'
import { NodeData } from 'const/NodeData'

export const ELEMENT_SLICE_NAME = 'element'

export type Param = {
  [name: string]: unknown
}

export type Algorithm = {
  [id: string]: {
    name: string
    param: Param
  }
}

export interface Element {
  flowElements: Elements<NodeData>
  clickedNodeId: string | null
  currentAlgoId: string
  algoParams: Algorithm
}
