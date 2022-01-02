import { RootState } from '../../store'

export const selectNwb = (state: RootState) => state.nwb

export const selectNwbList = (state: RootState) => {
  return selectNwb(state).nwbList
}
