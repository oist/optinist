import { RootState } from '../../store'

export const currentOutputIdSelector = (state: RootState) =>
  state.output.currentOutputId

export const currentOutputDataSelector = (state: RootState) => {
  const id = currentOutputIdSelector(state)
  if (Object.keys(state.output.outputData).includes(id)) {
    return state.output.outputData[id]
  } else {
    return undefined
  }
}
