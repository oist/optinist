import { enqueueSnackbar } from "notistack"

import { createSlice } from "@reduxjs/toolkit"

import { FILE_TREE_TYPE_SET } from "api/files/Files"
import { FileNodeFactory } from "factories/FileNodeFactory"
import { getFilesTree, deleteFile } from "store/slice/FilesTree/FilesTreeAction"
import {
  FilesTree,
  FILES_TREE_SLICE_NAME,
  TreeNodeType,
} from "store/slice/FilesTree/FilesTreeType"
import { convertToTreeNodeType } from "store/slice/FilesTree/FilesTreeUtils"
import { uploadFile } from "store/slice/FileUploader/FileUploaderActions"
// importSampleData available if needed for future workflow features

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
      .addCase(deleteFile.pending, (state, action) => {
        const { fileType } = action.meta.arg
        state[fileType] = {
          ...state[fileType],
          isLoading: true,
          isLatest: false,
        }
      })
      .addCase(deleteFile.rejected, (state, action) => {
        const { fileType } = action.meta.arg
        state[fileType] = {
          ...state[fileType],
          isLoading: false,
          isLatest: true,
        }
        enqueueSnackbar("Failed to delete file", { variant: "error" })
      })
      .addCase(deleteFile.fulfilled, (state, action) => {
        const { fileType, fileName } = action.meta.arg
        const filterTree = (nodes: TreeNodeType[]): TreeNodeType[] =>
          nodes
            .map((node) =>
              node.isDir ? { ...node, nodes: filterTree(node.nodes) } : node,
            )
            .filter((node) => node.path !== fileName)

        state[fileType].tree = filterTree(state[fileType].tree)
        state[fileType].isLoading = false
        state[fileType].isLatest = true
      })
      .addCase(uploadFile.pending, (state, action) => {
        const { fileType } = action.meta.arg
        // Get tree type from FileNodeFactory, fallback to ALL
        let treeType: string
        try {
          treeType = fileType
            ? FileNodeFactory.getTreeType(fileType)
            : FILE_TREE_TYPE_SET.ALL
        } catch (error) {
          // Fallback to ALL type for unknown file types
          treeType = FILE_TREE_TYPE_SET.ALL
        }

        state[treeType] = {
          ...state[treeType],
          isLoading: true,
          isLatest: false,
        }
      })
      .addCase(uploadFile.fulfilled, (state, action) => {
        const { fileType } = action.meta.arg

        // Use FileNodeFactory to get tree type dynamically
        if (fileType) {
          try {
            const config = FileNodeFactory.getFileTypeConfig(fileType)
            if (config && config.treeType && state[config.treeType]) {
              state[config.treeType].isLatest = false
            } else {
              state[FILE_TREE_TYPE_SET.ALL].isLatest = false
            }
          } catch {
            state[FILE_TREE_TYPE_SET.ALL].isLatest = false
          }
        } else {
          state[FILE_TREE_TYPE_SET.ALL].isLatest = false
        }
      })
  },
})

export default filesTreeSlice.reducer
