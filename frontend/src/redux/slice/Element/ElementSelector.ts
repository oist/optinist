import { RootState } from '../../store'

export const flowElementsSelector = (state: RootState) =>
  state.element.flowElements

export const currentElementSelector = (state: RootState) =>
  state.element.currentElement

export const algoParamsSelector = (state: RootState) => state.element.algoParams

export const paramValueSelector =
  (elementName: string, paramName: string) => (state: RootState) =>
    algoParamsSelector(state)[elementName][paramName]
