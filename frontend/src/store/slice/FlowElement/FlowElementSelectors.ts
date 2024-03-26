import { RootState } from "store/store"

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

const selectParentNodeIdsById = (nodeId: string) => (state: RootState) =>
  selectFlowEdges(state)
    .filter((edge) => edge.target === nodeId)
    .map((edge) => edge.source)

const selectNodeIsUpdate = (nodeId: string) => (state: RootState) =>
  state.algorithmNode[nodeId]?.isUpdate

export const selectAncestorNodesOriginalValueById =
  (nodeId: string) =>
  (state: RootState): boolean => {
    const parentNodes = selectParentNodeIdsById(nodeId)(state)
    if (parentNodes.length === 0) {
      return false
    } else {
      if (parentNodes.some((id) => selectNodeIsUpdate(id)(state))) {
        return true
      } else {
        return parentNodes.some((id) =>
          selectAncestorNodesOriginalValueById(id)(state),
        )
      }
    }
  }
