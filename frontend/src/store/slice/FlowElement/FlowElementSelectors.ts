import { RootState } from 'store/store'

export const selectFlowNodes = (state: RootState) => state.flowElement.flowNodes

export const selectFlowEdges = (state: RootState) => state.flowElement.flowEdges

export const selectFlowPosition = (state: RootState) =>
  state.flowElement.flowPosition

export const selectElementCoord = (state: RootState) =>
  state.flowElement.elementCoord

export const selectNodeById = (nodeId: string) => (state: RootState) =>
  selectFlowNodes(state).find((node) => node.id === nodeId)

export const selectNodeTypeById = (nodeId: string) => (state: RootState) =>
  selectNodeById(nodeId)(state)?.data?.type

export const selectNodeLabelById = (nodeId: string) => (state: RootState) =>
  selectNodeById(nodeId)(state)?.data?.label
