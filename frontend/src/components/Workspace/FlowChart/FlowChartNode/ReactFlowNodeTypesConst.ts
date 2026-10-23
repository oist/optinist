import React from "react"
import { NodeProps } from "reactflow"

import { CustomEdge } from "components/Workspace/FlowChart/CustomEdge"
import { getNodeComponent } from "components/Workspace/FlowChart/FlowChartNode/NodeComponentRegistry"
import { getAllFileTypeConfigs } from "config/fileTypes.config"
import { NodeData } from "store/slice/FlowElement/FlowElementType"

type ComponentType = React.ComponentType<NodeProps<NodeData>>

// Create node types mapping from config
const createNodeTypesFromConfig = () => {
  const nodeTypes: Record<string, ComponentType> = {}

  // Add AlgorithmNode (not part of file types)
  const algorithmNode = getNodeComponent("AlgorithmNode")
  if (algorithmNode) {
    nodeTypes.AlgorithmNode = algorithmNode
  }

  // Dynamically add file node types from config
  getAllFileTypeConfigs().forEach((config) => {
    const component = getNodeComponent(config.nodeComponent)
    if (component) {
      nodeTypes[config.reactFlowNodeType] = component
    }
  })

  return nodeTypes
}

export const reactFlowNodeTypes = createNodeTypesFromConfig()

export const reactFlowEdgeTypes = {
  buttonedge: CustomEdge,
} as const
