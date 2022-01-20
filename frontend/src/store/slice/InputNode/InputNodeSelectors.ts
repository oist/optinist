import { RootState } from 'store/store'

export const selectInputNode = (state: RootState) => state.inputNode

export const selectInputNodeById = (nodeId: string) => (state: RootState) =>
  state.inputNode[nodeId]

export const selectInputNodeDefined = (nodeId: string) => (state: RootState) =>
  Object.keys(state.inputNode).includes(nodeId)
export const selectInputNodeFileType = (nodeId: string) => (state: RootState) =>
  selectInputNodeById(nodeId)(state).fileType

export const selectInputNodeSelectedFilePath =
  (nodeId: string) => (state: RootState) =>
    selectInputNodeById(nodeId)(state).selectedFilePath

export const selectFilePathIsUndefined = (state: RootState) =>
  Object.values(state.inputNode).filter(
    (inputNode) => inputNode.selectedFilePath === undefined,
  ).length > 0
