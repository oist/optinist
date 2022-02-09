import { createAsyncThunk } from '@reduxjs/toolkit'
import axios from 'axios'
import { BASE_URL } from 'const/API'

import {
  FILES_TREE_SLICE_NAME,
  FILE_TREE_TYPE,
  TreeNodeTypeDTO,
} from './hdf5TreeType'

export const getHDF5FileTree = createAsyncThunk<
  TreeNodeTypeDTO[],
  FILE_TREE_TYPE
>(`${FILES_TREE_SLICE_NAME}/getFilesTree`, async (fileType, thunkAPI) => {
  try {
    const response = await axios.get(`${BASE_URL}/hdf5`, {
      params: {
        file_type: fileType,
      },
    })
    return response.data
  } catch (e) {
    return thunkAPI.rejectWithValue(e)
  }
})
