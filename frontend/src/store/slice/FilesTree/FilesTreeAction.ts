import { createAsyncThunk } from '@reduxjs/toolkit'
import axios from 'axios'
import { BASE_URL } from 'const/API'

import {
  FILES_TREE_SLICE_NAME,
  FILE_TYPE,
  TreeNodeTypeDTO,
} from './FilesTreeType'

export const getFilesTree = createAsyncThunk<TreeNodeTypeDTO[], FILE_TYPE>(
  `${FILES_TREE_SLICE_NAME}/getFilesTree`,
  async (fileType, thunkAPI) => {
    try {
      const response = await axios.get(`${BASE_URL}/files`, {
        params: {
          file_type: fileType,
        },
      })
      return response.data
    } catch (e) {
      return thunkAPI.rejectWithValue(e)
    }
  },
)
