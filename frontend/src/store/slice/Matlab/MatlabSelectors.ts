import { RootState } from "store/store"

export const selectMatlab = (state: RootState) => {
  if (state.hdf5 != null) {
    return state.hdf5
  } else {
    return undefined
  }
}

export const selectMatlabNodes = () => (state: RootState) =>
  selectMatlab(state)?.tree

export const selectMatlabIsLoading = () => (state: RootState) =>
  selectMatlab(state)?.isLoading ?? false
