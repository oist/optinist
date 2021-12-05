import { AlgoNodeData } from 'const/NodeData'
import { ElementId } from 'react-flow-renderer'
import { isNodeData, isInputNodeData, isAlgoNodeData } from 'utils/ElementUtils'
import { RootState } from '../../store'
import { algoParamByIdSelector } from '../Algorithm/AlgorithmSelector'

export const flowElementsSelector = (state: RootState) =>
  state.element.flowElements

export const nodeByIdSelector =
  (elementId: ElementId) => (state: RootState) => {
    return flowElementsSelector(state)
      .filter(isNodeData)
      .find((node) => node.id === elementId)
  }

export const maxElementIdSelector = (state: RootState) =>
  state.element.flowElements
    .map((element) => element.id)
    .map((id) => Number(id))
    .filter((id) => !isNaN(id))
    .reduce((a, b) => Math.max(a, b))

export const nodeDataListForRunSelector = (state: RootState) =>
  flowElementsSelector(state)
    .filter((element) => isInputNodeData(element) || isAlgoNodeData(element))
    .map((element) => {
      if (element.data && element.data.type === 'algo') {
        const param = algoParamByIdSelector(element.id)(state)
        const data: AlgoNodeData = {
          ...element.data,
          param,
        }
        return data
      } else {
        return element.data
      }
    })

export const pathIsUndefinedSelector = (state: RootState) => {
  const pathErrorNodeList = flowElementsSelector(state).filter((element) => {
    if (isInputNodeData(element)) {
      if (!element.data?.path) {
        return true
      }
    }
    return false
  })
  return pathErrorNodeList.length > 0
}

export const filePathSelector = (nodeId: string) => (state: RootState) => {
  const node = nodeByIdSelector(nodeId)(state)
  return node?.data?.path
}
