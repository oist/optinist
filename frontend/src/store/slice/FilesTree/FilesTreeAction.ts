import { createAsyncThunk } from '@reduxjs/toolkit'
import {
  FILE_TREE_TYPE,
  getFilesTreeApi,
  TreeNodeTypeDTO,
} from 'api/files/Files'

import { FILES_TREE_SLICE_NAME } from './FilesTreeType'

export const getFilesTree = createAsyncThunk<TreeNodeTypeDTO[], FILE_TREE_TYPE>(
  `${FILES_TREE_SLICE_NAME}/getFilesTree`,
  async (fileType, thunkAPI) => {
    try {
      const response = await getFilesTreeApi(fileType)
      return response
    } catch (e) {
      return thunkAPI.rejectWithValue(e)
    }
  },
)
