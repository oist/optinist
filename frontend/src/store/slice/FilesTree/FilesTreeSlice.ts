import { createSlice } from '@reduxjs/toolkit'
import { uploadFile } from 'store/slice/FileUploader/FileUploaderActions'
import { FILE_TYPE_SET } from 'store/slice/InputNode/InputNodeType'

import { getFilesTree } from './FilesTreeAction'
import {
  FilesTree,
  FILES_TREE_SLICE_NAME,
  FILE_TREE_TYPE_SET,
} from './FilesTreeType'
import { convertToTreeNodeType } from './FilesTreeUtils'

const initialState: FilesTree = {}
export const filesTreeSlice = createSlice({
  name: FILES_TREE_SLICE_NAME,
  initialState,
  reducers: {},
  extraReducers: (builder) => {
    builder
      .addCase(getFilesTree.pending, (state, action) => {
        const fileType = action.meta.arg
        state[fileType] = {
          isLoading: true,
          isLatest: false,
          tree: [],
        }
      })
      .addCase(getFilesTree.fulfilled, (state, action) => {
        const fileTreeType = action.meta.arg
        state[fileTreeType].tree = convertToTreeNodeType(action.payload)
        state[fileTreeType].isLatest = true
        state[fileTreeType].isLoading = false
      })
      .addCase(uploadFile.pending, (state, action) => {
        const { fileType } = action.meta.arg
        if (fileType === FILE_TYPE_SET.IMAGE) {
          if (state[FILE_TREE_TYPE_SET.IMAGE] != null) {
            state[FILE_TREE_TYPE_SET.IMAGE].isLatest = false
          } else {
            state[FILE_TREE_TYPE_SET.IMAGE] = {
              isLoading: false,
              isLatest: false,
              tree: [],
            }
          }
        } else if (fileType === FILE_TYPE_SET.CSV) {
          if (state[FILE_TREE_TYPE_SET.CSV] != null) {
            state[FILE_TREE_TYPE_SET.CSV].isLatest = false
          } else {
            state[FILE_TREE_TYPE_SET.CSV] = {
              isLoading: false,
              isLatest: false,
              tree: [],
            }
          }
        } else {
          if (state[FILE_TREE_TYPE_SET.ALL] != null) {
            state[FILE_TREE_TYPE_SET.ALL].isLatest = false
          } else {
            state[FILE_TREE_TYPE_SET.ALL] = {
              isLoading: false,
              isLatest: false,
              tree: [],
            }
          }
        }
      })
      .addCase(uploadFile.fulfilled, (state, action) => {
        const { fileType } = action.meta.arg
        if (fileType === FILE_TYPE_SET.IMAGE) {
          state[FILE_TREE_TYPE_SET.IMAGE].isLatest = false
        } else if (fileType === FILE_TYPE_SET.CSV) {
          state[FILE_TREE_TYPE_SET.CSV].isLatest = false
        } else {
          state[FILE_TREE_TYPE_SET.ALL].isLatest = false
        }
      })
  },
})

export default filesTreeSlice.reducer
