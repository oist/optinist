import { createSlice } from '@reduxjs/toolkit'

import { uploadCsvFile, uploadImageFile } from '../FileData/FileDataAction'
import { getFilesTree } from './FilesTreeAction'
import {
  FilesTree,
  FILES_TREE_SLICE_NAME,
  FILE_TYPE_SET,
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
        const fileType = action.meta.arg
        state[fileType].tree = convertToTreeNodeType(action.payload)
        state[fileType].isLatest = true
        state[fileType].isLoading = false
      })
      .addCase(uploadImageFile.pending, (state) => {
        if (state[FILE_TYPE_SET.IMAGE] != null) {
          state[FILE_TYPE_SET.IMAGE].isLatest = false
        } else {
          state[FILE_TYPE_SET.IMAGE] = {
            isLoading: false,
            isLatest: false,
            tree: [],
          }
        }
      })
      .addCase(uploadImageFile.fulfilled, (state) => {
        state[FILE_TYPE_SET.IMAGE].isLatest = false
      })
      .addCase(uploadCsvFile.pending, (state) => {
        if (state[FILE_TYPE_SET.CSV] != null) {
          state[FILE_TYPE_SET.CSV].isLatest = false
        } else {
          state[FILE_TYPE_SET.CSV] = {
            isLoading: false,
            isLatest: false,
            tree: [],
          }
        }
      })
      .addCase(uploadCsvFile.fulfilled, (state) => {
        state[FILE_TYPE_SET.CSV].isLatest = false
      })
  },
})

export default filesTreeSlice.reducer
