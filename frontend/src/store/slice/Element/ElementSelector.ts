import { ElementId } from 'react-flow-renderer'
import { isNodeData, isImageNodeData } from 'utils/ElementUtils'
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
  state.element.flowElements.length === 0
    ? 0
    : state.element.flowElements
        .map((element) => element.id)
        .map((id) => Number(id))
        .filter((id) => !isNaN(id))
        .reduce((a, b) => Math.max(a, b))

export const nodeDataListForRunSelector = (state: RootState) =>
  flowElementsSelector(state).map((element) => {
    if (element.data && element.data.type === 'algo') {
      const param = algoParamByIdSelector(element.id)(state)
      return {
        ...element,
        data: {
          ...element.data,
          param,
        },
      }
    } else {
      return element
    }
  })

export const pathIsUndefinedSelector = (state: RootState) => {
  const pathErrorNodeList = flowElementsSelector(state).filter((element) => {
    if (isImageNodeData(element)) {
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
