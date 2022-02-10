import { RootState } from '../../store'

export const selectHDF5 = (state: RootState) => {
  if (state.filesTree != null) {
    return state.filesTree
  } else {
    return undefined
  }
}

export const selectHDF5Nodes = () => (state: RootState) =>
  selectHDF5(state)?.tree

export const selectHDF5IsLatest = () => (state: RootState) =>
  selectHDF5(state)?.isLatest ?? false

export const selectHDF5IsLoading = () => (state: RootState) =>
  selectHDF5(state)?.isLoading ?? false
