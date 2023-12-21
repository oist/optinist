import { RootState } from "store/store"

export const selectMatlab = (state: RootState) => {
  if (state.matlab != null) {
    return state.matlab
  } else {
    return undefined
  }
}

export const selectMatlabNodes = () => (state: RootState) =>
  selectMatlab(state)?.tree

export const selectMatlabIsLoading = () => (state: RootState) =>
  selectMatlab(state)?.isLoading ?? false
