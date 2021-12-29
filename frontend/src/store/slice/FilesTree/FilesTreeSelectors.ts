import { RootState } from '../../store'
import { FILE_TREE_TYPE } from './FilesTreeType'

export const selectFilesTree =
  (fileType: FILE_TREE_TYPE) => (state: RootState) => {
    if (state.filesTree[fileType] != null) {
      return state.filesTree[fileType]
    } else {
      return undefined
    }
  }

export const selectFilesTreeNodes =
  (fileType: FILE_TREE_TYPE) => (state: RootState) =>
    selectFilesTree(fileType)(state)?.tree

export const selectFilesIsLatest =
  (fileType: FILE_TREE_TYPE) => (state: RootState) =>
    selectFilesTree(fileType)(state)?.isLatest ?? false

export const selectFilesIsLoading =
  (fileType: FILE_TREE_TYPE) => (state: RootState) =>
    selectFilesTree(fileType)(state)?.isLoading ?? false
