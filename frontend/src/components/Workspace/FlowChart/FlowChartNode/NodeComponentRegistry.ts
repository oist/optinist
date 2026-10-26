import React from "react"
import { NodeProps } from "reactflow"

// Import all node components
import { AlgorithmNode } from "components/Workspace/FlowChart/FlowChartNode/AlgorithmNode"
import { BehaviorFileNode } from "components/Workspace/FlowChart/FlowChartNode/BehaviorFileNode"
import { CsvFileNode } from "components/Workspace/FlowChart/FlowChartNode/CsvFileNode"
import { FluoFileNode } from "components/Workspace/FlowChart/FlowChartNode/FluoFileNode"
import { HDF5FileNode } from "components/Workspace/FlowChart/FlowChartNode/HDF5FileNode"
import { ImageFileNode } from "components/Workspace/FlowChart/FlowChartNode/ImageFileNode"
import { MatlabFileNode } from "components/Workspace/FlowChart/FlowChartNode/MatlabFileNode"
import { MicroscopeFileNode } from "components/Workspace/FlowChart/FlowChartNode/MicroscopeFileNode"
import { NodeData } from "store/slice/FlowElement/FlowElementType"

type NodeComponentType = React.ComponentType<NodeProps<NodeData>>

// Component registry mapping node type names to components
export const nodeComponentRegistry: Record<string, NodeComponentType> = {
  AlgorithmNode,
  ImageFileNode,
  CsvFileNode,
  MatlabFileNode,
  HDF5FileNode,
  FluoFileNode,
  BehaviorFileNode,
  MicroscopeFileNode,
}

// Get component by node type name
export function getNodeComponent(
  nodeType: string,
): NodeComponentType | undefined {
  return nodeComponentRegistry[nodeType]
}

// Check if component exists
export function hasNodeComponent(nodeType: string): boolean {
  return nodeType in nodeComponentRegistry
}
