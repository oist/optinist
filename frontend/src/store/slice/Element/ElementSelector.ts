import { ElementId } from 'react-flow-renderer'
import { isNodeData } from 'utils/ElementUtils'
import { RootState } from '../../store'

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

export const runStatusSelector = (state: RootState) => state.element.runStatus

export const runMassageSelector = (state: RootState) => state.element.runMessage
