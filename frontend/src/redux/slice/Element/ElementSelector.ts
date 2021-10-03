import { ElementId } from 'react-flow-renderer'
import { RootState } from '../../store'

export const flowElementsSelector = (state: RootState) =>
  state.element.flowElements

export const currentAlgoIdSelector = (state: RootState) =>
  state.element.currentAlgoId

export const elementByIdSelector =
  (elementId: ElementId) => (state: RootState) => {
    return state.element.flowElements.find((node) => node.id === elementId)
  }

export const algoParamsSelector = (state: RootState) => state.element.algoParams

export const algoParamByIdSelector = (id: string) => (state: RootState) => {
  const algoParams = algoParamsSelector(state)
  if (Object.keys(algoParams).includes(id)) {
    return algoParams[id]
  } else {
    return undefined
  }
}

export const paramValueSelector =
  (elementId: string, paramName: string) => (state: RootState) =>
    algoParamsSelector(state)[elementId].param[paramName]
