import { Elements } from 'react-flow-renderer'
import { NodeData, NODE_DATA_TYPE_SET } from 'const/NodeData'

export const INITIAL_IMAGE_ELEMENT_ID = '0'
export const INITIAL_IMAGE_ELEMENT_NAME = 'ImageData'
export const INITIAL_ALGO_ELEMENT_ID = '1'
export const INITIAL_ALGO_ELEMENT_NAME = 'dummy_image2image'
export const INITIAL_EDGE_ID = '2'

export const INITIAL_ALGO_STYLE: React.CSSProperties = {
  width: 180,
  height: 100,
  padding: 0,
  borderRadius: 0,
} as const

export const INITIAL_DATA_STYLE: React.CSSProperties = {
  border: '1px solid #777',
  height: 120,
} as const

export const initialElements: Elements<NodeData> = [
  {
    id: INITIAL_IMAGE_ELEMENT_ID,
    type: 'ImageFileNode',
    data: {
      type: NODE_DATA_TYPE_SET.IMAGE,
      label: INITIAL_IMAGE_ELEMENT_NAME,
    },
    style: INITIAL_DATA_STYLE,
    position: { x: 50, y: 60 },
  },
  // {
  //   id: INITIAL_ALGO_ELEMENT_ID,
  //   type: 'default',
  //   data: { type: NODE_DATA_TYPE_SET.ALGO, label: INITIAL_ALGO_ELEMENT_NAME },
  //   style: INITIAL_ALGO_STYLE,
  //   position: { x: 400, y: 60 },
  //   targetPosition: Position.Left,
  //   sourcePosition: Position.Right,
  // },
  // // edge
  // {
  //   id: INITIAL_EDGE_ID,
  //   source: INITIAL_IMAGE_ELEMENT_ID,
  //   target: INITIAL_ALGO_ELEMENT_ID,
  //   targetHandle: `${INITIAL_ALGO_ELEMENT_ID}-image`,
  //   type: 'smoothstep',
  // },
]
