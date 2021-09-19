const initialElements = [
  // {
  //   id: '1',
  //   type: 'input',
  //   data: { label: 'data1' },
  //   position: { x: 200, y: 5 },
  // },
  {
    id: '1',
    type: 'selectorNode',
    data: { label: 'data' },
    style: { border: '1px solid #777', padding: 10 },
    position: { x: 200, y: 5 },
  },
  {
    id: '2',
    type: 'default',
    data: { label: 'caiman_mc' },
    position: { x: 200, y: 100 },
  },
  {
    id: '3',
    type: 'default',
    data: { label: 'caiman_cnmf' },
    position: { x: 200, y: 200 },
  },
  {
    id: '4',
    type: 'output',
    data: { label: 'output' },
    position: { x: 200, y: 300 },
  },
  // edge
  {
    id: '1',
    source: '1',
    target: '2',
    type: 'smoothstep',
  },
  {
    id: '2',
    source: '2',
    target: '3',
    type: 'smoothstep',
  },
  {
    id: '3',
    source: '3',
    target: '4',
    type: 'smoothstep',
  },
]

export { initialElements }
