import { RootState } from '../../store'
import { FILE_TYPE } from './FilesTreeType'

export const filesTreeSelector =
  (fileType: FILE_TYPE) => (state: RootState) => {
    if (state.filesTree[fileType] != null) {
      return state.filesTree[fileType]
    } else {
      return undefined
    }
  }

export const filesTreeNodesSelector =
  (fileType: FILE_TYPE) => (state: RootState) =>
    filesTreeSelector(fileType)(state)?.tree

export const filesIsLatestSelector =
  (fileType: FILE_TYPE) => (state: RootState) =>
    filesTreeSelector(fileType)(state)?.isLatest ?? false

export const filesIsLoadingSelector =
  (fileType: FILE_TYPE) => (state: RootState) =>
    filesTreeSelector(fileType)(state)?.isLoading ?? false
