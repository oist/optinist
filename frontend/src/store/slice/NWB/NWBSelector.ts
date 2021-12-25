import { RootState } from '../../store'

export const nwbSelector = (state: RootState) => state.nwb

export const nwbListSelector = (state: RootState) => {
  return nwbSelector(state).nwbList
}
