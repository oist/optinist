import { Elements } from 'react-flow-renderer'
import { NodeData } from 'const/NodeData'

export const ELEMENT_SLICE_NAME = 'element'

export interface Element {
  flowElements: Elements<NodeData>
}
