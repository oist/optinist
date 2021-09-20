import { RootState } from '../../store'

export const flowElementsSelector = (state: RootState) =>
  state.element.flowElements
