import { isEdge } from 'react-flow-renderer'
import { RootState } from 'store/store'
import {
  selectAlgorithmFunctionPath,
  selectAlgorithmParams,
} from '../AlgorithmNode/AlgorithmNodeSelectors'
import {
  selectInputNodeHDF5Path,
  selectInputNodeSelectedFilePath,
} from '../InputNode/InputNodeSelectors'
import { NODE_TYPE_SET } from './FlowElementType'
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

export const selectElementListForRun = (state: RootState) => {
  const elements = selectFlowElements(state)
  const nodeList = elements.filter(isNodeData).map((element) => {
    if (element.data) {
      if (element.data.type === NODE_TYPE_SET.ALGORITHM) {
        const param = selectAlgorithmParams(element.id)(state) ?? {}
        const functionPath = selectAlgorithmFunctionPath(element.id)(state)
        return {
          ...element,
          data: {
            ...element.data,
            param,
            path: functionPath,
          },
        }
      } else if (element.data.type === NODE_TYPE_SET.INPUT) {
        const filePath = selectInputNodeSelectedFilePath(element.id)(state)
        const hdf5Path = selectInputNodeHDF5Path(element.id)(state)
        return {
          ...element,
          data: {
            ...element.data,
            path: filePath,
            hdf5Path: hdf5Path,
          },
        }
      }
    }
    return element
  })
  const edgeList = elements.filter(isEdge)
  return { nodeList, edgeList }
}
