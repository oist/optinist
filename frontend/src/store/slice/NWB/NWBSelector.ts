import { RootState } from '../../store'

export const nwbSelector = (state: RootState) => state.nwb

export const nwbParamValueSelector =
  (paramName: string) => (state: RootState) => {
    const value = nwbSelector(state).params[paramName]
    return value
  }

export const nwbParamNameListSelector = (state: RootState) => {
  return Object.keys(nwbSelector(state).params)
}

export const nwbListSelector = (state: RootState) => {
  return nwbSelector(state).nwbList
}
