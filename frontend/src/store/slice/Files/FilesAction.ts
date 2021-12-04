import { createAsyncThunk } from '@reduxjs/toolkit'
import axios from 'axios'
import { BASE_URL } from 'const/API'

import { FILES_SLICE_NAME, TreeNodeTypeDTO } from './FilesType'

export const getFiles = createAsyncThunk<TreeNodeTypeDTO[], string | undefined>(
  `${FILES_SLICE_NAME}/getFiles`,
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
