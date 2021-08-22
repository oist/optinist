const initialElements = [
  {
    id: '1',
    type: 'input',
    data: { label: 'data1' },
    position: { x: 200, y: 5 },
  },
  {
    id: '2',
    type: 'default',
    data: { label: 'algo1' },
    position: { x: 200, y: 100 },
  },
  {
    id: '3',
    type: 'output',
    data: { label: 'output1' },
    position: { x: 200, y: 200 },
  },
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
]

export { initialElements }
