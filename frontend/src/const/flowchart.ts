import { Elements } from 'react-flow-renderer'
import { NodeData, NODE_DATA_TYPE_SET } from 'const/NodeData'

export const INITIAL_IMAGE_ELEMENT_ID = '0'
export const INITIAL_ALGO_ELEMENT_ID = '1'
export const INITIAL_OUTPUT_ELEMENT_ID = '2'

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
    style: { border: '1px solid #777', padding: 10 },
    position: { x: 200, y: 5 },
  },
  {
    id: INITIAL_ALGO_ELEMENT_ID,
    type: 'default',
    data: { type: NODE_DATA_TYPE_SET.ALGO, label: 'caiman_mc' },
    position: { x: 200, y: 100 },
  },
  // {
  //   id: '3',
  //   type: 'default',
  //   data: { type: 'algo', label: 'caiman_cnmf' },
  //   position: { x: 200, y: 200 },
  // },
  {
    id: INITIAL_OUTPUT_ELEMENT_ID,
    type: 'output',
    data: { type: NODE_DATA_TYPE_SET.OUTPUT, label: 'output' },
    position: { x: 200, y: 300 },
  },

  // edge
  {
    id: '3',
    source: INITIAL_IMAGE_ELEMENT_ID,
    target: INITIAL_ALGO_ELEMENT_ID,
    type: 'smoothstep',
  },
  {
    id: '4',
    source: INITIAL_ALGO_ELEMENT_ID,
    target: INITIAL_OUTPUT_ELEMENT_ID,
    type: 'smoothstep',
  },
  // {
  //   id: 'e3',
  //   source: '3',
  //   target: '4',
  //   type: 'smoothstep',
  // },
]
