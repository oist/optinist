import { createAsyncThunk } from "@reduxjs/toolkit"

import {
  FILE_TREE_TYPE,
  getFilesTreeApi,
  deleteFileApi,
  TreeNodeTypeDTO,
} from "api/files/Files"
import { FILES_TREE_SLICE_NAME } from "store/slice/FilesTree/FilesTreeType"

export const getFilesTree = createAsyncThunk<
  TreeNodeTypeDTO[],
  { workspaceId: number; fileType: FILE_TREE_TYPE }
>(
  `${FILES_TREE_SLICE_NAME}/getFilesTree`,
  async ({ workspaceId, fileType }, thunkAPI) => {
    try {
      const response = await getFilesTreeApi(workspaceId, fileType)
      return response
    } catch (e) {
      return thunkAPI.rejectWithValue(e)
    }
  },
)

export const deleteFile = createAsyncThunk<
  boolean,
  { workspaceId: number; fileName: string; fileType: FILE_TREE_TYPE }
>(
  `${FILES_TREE_SLICE_NAME}/deleteFile`,
  async ({ workspaceId, fileName }, thunkAPI) => {
    if (workspaceId) {
      try {
        const response = await deleteFileApi(workspaceId, fileName)
        return response
      } catch (e) {
        return thunkAPI.rejectWithValue(e)
      }
    } else {
      return thunkAPI.rejectWithValue("workspace id does not exist.")
    }
  },
)
