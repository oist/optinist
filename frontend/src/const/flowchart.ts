import { Elements } from 'react-flow-renderer'

export const INITIAL_ALGO_ELEMENT_ID = '2'

export const initialElements: Elements = [
  // {
  //   id: '1',
  //   type: 'input',
  //   data: { label: 'data1' },
  //   position: { x: 200, y: 5 },
  // },
  {
    id: '1',
    type: 'selectorNode',
    data: {
      type: 'data',
      label: 'data',
      path: '/Users/shogoakiyama/caiman_data/example_movies/Sue_2x_3000_40_-46.tif',
    },
    style: { border: '1px solid #777', padding: 10 },
    position: { x: 200, y: 5 },
  },
  {
    id: INITIAL_ALGO_ELEMENT_ID,
    type: 'default',
    data: { type: 'algo', label: 'caiman_mc' },
    position: { x: 200, y: 100 },
  },
  // {
  //   id: '3',
  //   type: 'default',
  //   data: { type: 'algo', label: 'caiman_cnmf' },
  //   position: { x: 200, y: 200 },
  // },
  {
    id: '3',
    type: 'output',
    data: { type: 'output', label: 'output' },
    position: { x: 200, y: 300 },
  },

  // edge
  {
    id: 'e1',
    source: '1',
    target: INITIAL_ALGO_ELEMENT_ID,
    type: 'smoothstep',
  },
  {
    id: 'e2',
    source: INITIAL_ALGO_ELEMENT_ID,
    target: '3',
    type: 'smoothstep',
  },
  // {
  //   id: 'e3',
  //   source: '3',
  //   target: '4',
  //   type: 'smoothstep',
  // },
]
