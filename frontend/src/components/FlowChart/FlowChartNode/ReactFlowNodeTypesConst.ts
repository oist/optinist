import { ImageFileNode } from 'components/FlowChart/FlowChartNode/ImageFileNode'
import { AlgorithmNode } from 'components/FlowChart/FlowChartNode/AlgorithmNode'
import { CsvFileNode } from 'components/FlowChart/FlowChartNode/CsvFileNode'
import { HDF5FileNode } from 'components/FlowChart/FlowChartNode/HDF5FileNode'
import { FluoFileNode } from 'components/FlowChart/FlowChartNode/FluoFileNode'
import { BehaviorFileNode } from 'components/FlowChart/FlowChartNode/BehaviorFileNode'

import { CustomEdge } from 'components/FlowChart/CustomEdge'

export const reactFlowNodeTypes = {
  ImageFileNode,
  CsvFileNode,
  HDF5FileNode,
  AlgorithmNode,
  FluoFileNode,
  BehaviorFileNode,
} as const

export const reactFlowEdgeTypes = {
  buttonedge: CustomEdge,
} as const
