import { createSlice } from "@reduxjs/toolkit"

import { FILE_TREE_TYPE_SET } from "api/files/Files"
import {
  getFilesTree,
  deleteFileTree,
} from "store/slice/FilesTree/FilesTreeAction"
import {
  FilesTree,
  FILES_TREE_SLICE_NAME,
} from "store/slice/FilesTree/FilesTreeType"
import { convertToTreeNodeType } from "store/slice/FilesTree/FilesTreeUtils"
import { uploadFile } from "store/slice/FileUploader/FileUploaderActions"
import { FILE_TYPE_SET } from "store/slice/InputNode/InputNodeType"
import { importSampleData } from "store/slice/Workflow/WorkflowActions"

export const initialState: FilesTree = {}
export const filesTreeSlice = createSlice({
  name: FILES_TREE_SLICE_NAME,
  initialState,
  reducers: {},
  extraReducers: (builder) => {
    builder
      .addCase(getFilesTree.pending, (state, action) => {
        const { fileType } = action.meta.arg
        state[fileType] = {
          ...state[fileType],
          isLoading: true,
          isLatest: false,
        }
      })
      .addCase(getFilesTree.fulfilled, (state, action) => {
        const { fileType } = action.meta.arg
        state[fileType].tree = convertToTreeNodeType(action.payload)
        state[fileType].isLatest = true
        state[fileType].isLoading = false
      })
      .addCase(deleteFileTree.pending, (state, action) => {
        const { fileType } = action.meta.arg
        console.log(state[fileType].tree)
        console.log(state[fileType].isLatest)
        console.log(state[fileType].isLoading)
        console.log(action.meta.arg)
        state[fileType] = {
          ...state[fileType],
          isLoading: true,
          isLatest: false,
        }
      })
      .addCase(deleteFileTree.rejected, (state, action) => {
        const { fileType } = action.meta.arg
        state[fileType] = {
          ...state[fileType],
          isLoading: false,
          isLatest: true,
        }
        console.log("deleteFileTree.rejected")
        console.log(action.meta.arg)
        console.log(state[fileType].tree)
        console.log(state[fileType].isLatest)
        console.log(state[fileType].isLoading)
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
        } else if (fileType === FILE_TYPE_SET.HDF5) {
          if (state[FILE_TREE_TYPE_SET.HDF5] != null) {
            state[FILE_TREE_TYPE_SET.HDF5].isLatest = false
          } else {
            state[FILE_TREE_TYPE_SET.HDF5] = {
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
        } else if (fileType === FILE_TYPE_SET.HDF5) {
          state[FILE_TREE_TYPE_SET.HDF5].isLatest = false
        } else {
          state[FILE_TREE_TYPE_SET.ALL].isLatest = false
        }
      })
      .addCase(importSampleData.fulfilled, (state) => {
        ;[
          FILE_TREE_TYPE_SET.IMAGE,
          FILE_TREE_TYPE_SET.CSV,
          FILE_TREE_TYPE_SET.HDF5,
          FILE_TREE_TYPE_SET.ALL,
        ].forEach((fileType) => {
          if (state[fileType] != null) {
            state[fileType].isLatest = false
          }
        })
      })
  },
})

export default filesTreeSlice.reducer
