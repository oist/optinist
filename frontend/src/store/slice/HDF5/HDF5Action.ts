import { createAsyncThunk } from '@reduxjs/toolkit'
import axios from 'axios'
import { BASE_URL } from 'const/API'

import { HDF5_SLICE_NAME, TreeNodeTypeDTO } from './HDF5Type'

export const getHDF5Tree = createAsyncThunk<TreeNodeTypeDTO[], void>(
  `${HDF5_SLICE_NAME}/getHDF5Tree`,
  async (thunkAPI) => {
    try {
      const response = await axios.get(`${BASE_URL}/hdf5`)
      return response.data
    } catch (e) {
      // return thunkAPI.rejectWithValue(e)
    }
  },
)
