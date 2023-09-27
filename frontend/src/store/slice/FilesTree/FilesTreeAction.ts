import { createAsyncThunk } from '@reduxjs/toolkit'
import {
  FILE_TREE_TYPE,
  getFilesTreeApi,
  TreeNodeTypeDTO,
} from 'api/files/Files'

import { FILES_TREE_SLICE_NAME } from './FilesTreeType'

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
