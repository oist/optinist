import { Elements, Position } from 'react-flow-renderer'
import { NodeData, NODE_DATA_TYPE_SET } from 'const/NodeData'

export const INITIAL_IMAGE_ELEMENT_ID = '0'
export const INITIAL_ALGO_ELEMENT_ID = '1'
export const INITIAL_EDGE_ID = '2'
export const INITIAL_OUTPUT_ELEMENT_ID = '3'

export const INITIAL_ALGO_STYLE: React.CSSProperties = {
  width: 180,
  height: 100,
  padding: 0,
  borderRadius: 0,
} as const

export const INITIAL_DATA_STYLE: React.CSSProperties = {
  border: '1px solid #777',
  height: 100,
} as const

export const initialElements: Elements<NodeData> = [
  // {
  //   id: '1',
  //   type: 'input',
  //   data: { label: 'data1' },
  //   position: { x: 200, y: 5 },
  // },
  {
    id: INITIAL_IMAGE_ELEMENT_ID,
    type: 'selectorNode',
    data: {
      type: NODE_DATA_TYPE_SET.DATA,
      label: 'data',
      path: '/Users/shogoakiyama/caiman_data/example_movies/Sue_2x_3000_40_-46.tif',
    },
    style: INITIAL_DATA_STYLE,
    position: { x: 50, y: 60 },
  },
  {
    id: INITIAL_ALGO_ELEMENT_ID,
    type: 'default',
    data: { type: NODE_DATA_TYPE_SET.ALGO, label: 'caiman_mc' },
    style: INITIAL_ALGO_STYLE,
    position: { x: 400, y: 60 },
    targetPosition: Position.Left,
    sourcePosition: Position.Right,
  },
  // {
  //   id: '3',
  //   type: 'default',
  //   data: { type: 'algo', label: 'caiman_cnmf' },
  //   position: { x: 200, y: 200 },
  // },
  // {
  //   id: INITIAL_OUTPUT_ELEMENT_ID,
  //   type: 'output',
  //   data: { type: NODE_DATA_TYPE_SET.OUTPUT, label: 'output' },
  //   position: { x: 200, y: 300 },
  // },

  // edge
  {
    id: INITIAL_EDGE_ID,
    source: INITIAL_IMAGE_ELEMENT_ID,
    target: INITIAL_ALGO_ELEMENT_ID,
    type: 'smoothstep',
  },
  // {
  //   id: '4',
  //   source: INITIAL_ALGO_ELEMENT_ID,
  //   target: INITIAL_OUTPUT_ELEMENT_ID,
  //   type: 'smoothstep',
  // },
  // {
  //   id: 'e3',
  //   source: '3',
  //   target: '4',
  //   type: 'smoothstep',
  // },
]
