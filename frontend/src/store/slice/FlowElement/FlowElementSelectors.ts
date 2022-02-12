import { isEdge } from 'react-flow-renderer'
import { RootState } from 'store/store'

import { isNodeData } from './FlowElementUtils'

export const selectFlowElements = (state: RootState) =>
  state.flowElement.flowElements

export const selectFlowPosition = (state: RootState) =>
  state.flowElement.flowPosition

export const selectElementCoord = (state: RootState) =>
  state.flowElement.elementCoord

export const selectNodeById = (nodeId: string) => (state: RootState) =>
  selectFlowElements(state)
    .filter(isNodeData)
    .find((node) => node.id === nodeId)

export const selectNodeTypeById = (nodeId: string) => (state: RootState) =>
  selectNodeById(nodeId)(state)?.data?.type

export const selectNodeLabelById = (nodeId: string) => (state: RootState) =>
  selectNodeById(nodeId)(state)?.data?.label

export const selectEdgeListForRun = (state: RootState) => {
  return selectFlowElements(state).filter(isEdge)
}
