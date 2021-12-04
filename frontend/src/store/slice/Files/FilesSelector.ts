import { RootState } from '../../store'

export const filesSelector = (state: RootState) => state.files

export const filesTreeSelector = (state: RootState) => filesSelector(state).tree

export const filesIsLatestSelector = (state: RootState) =>
  filesSelector(state).isLatest

export const filesIsLoadingSelector = (state: RootState) =>
  filesSelector(state).isLoading
