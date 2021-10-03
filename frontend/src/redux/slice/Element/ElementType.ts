import { Elements } from 'react-flow-renderer'

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
  flowElements: Elements
  currentElementId: string
  algoParams: Algorithm
}

export type NodeData = {
  label: string
  path?: string
}
