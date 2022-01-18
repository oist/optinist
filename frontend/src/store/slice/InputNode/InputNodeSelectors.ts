import { RootState } from 'store/store'
import { selectNodeTypeById } from '../FlowElement/FlowElementSelectors'
import { NODE_TYPE_SET } from '../FlowElement/FlowElementType'
import { isImageInputNode } from './InputNodeUtils'

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

export const selectImageInputNodeMaxIndex =
  (nodeId: string) => (state: RootState) => {
    const inputNode = selectInputNodeById(nodeId)(state)
    if (isImageInputNode(inputNode)) {
      return inputNode.maxIndex
    } else {
      throw new Error('invalid nodeId of ImageInputNode')
    }
  }

export const selectImageMaxIndexByNodeId =
  (nodeId: string) => (state: RootState) => {
    const nodeType = selectNodeTypeById(nodeId)(state)
    if (nodeType === NODE_TYPE_SET.INPUT) {
      return selectImageInputNodeMaxIndex(nodeId)(state)
    } else {
      return null
    }
  }

export const selectFilePathIsUndefined = (state: RootState) =>
  Object.values(state.inputNode).filter(
    (inputNode) => inputNode.selectedFilePath === undefined,
  ).length > 0
