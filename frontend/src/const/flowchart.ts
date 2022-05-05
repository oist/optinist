export const INITIAL_IMAGE_ELEMENT_ID = 'input_0'
export const INITIAL_IMAGE_ELEMENT_NAME = 'NoName'

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

export const REACT_FLOW_NODE_TYPE_KEY = {
  ImageFileNode: 'ImageFileNode',
  CsvFileNode: 'CsvFileNode',
  HDF5FileNode: 'HDF5FileNode',
  FluoFileNode: 'FluoFileNode',
  AlgorithmNode: 'AlgorithmNode',
  BehaviorFileNode: 'BehaviorFileNode',
} as const

export type REACT_FLOW_NODE_TYPE =
  typeof REACT_FLOW_NODE_TYPE_KEY[keyof typeof REACT_FLOW_NODE_TYPE_KEY]
