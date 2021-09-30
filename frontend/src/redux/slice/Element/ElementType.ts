import { Elements } from 'react-flow-renderer'

export type Param = {
  [name: string]: number
}

export type Algorithm = {
  [name: string]: Param
}

export interface Element {
  flowElements: Elements
  currentElement: string
  algoParams: Algorithm
}
