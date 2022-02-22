import { createAsyncThunk } from '@reduxjs/toolkit'
import axios from 'axios'
import { BASE_URL } from 'const/API'

import { HDF5_SLICE_NAME, HDF5TreeDTO } from './HDF5Type'

export const getHDF5Tree = createAsyncThunk<HDF5TreeDTO[], { path: string }>(
  `${HDF5_SLICE_NAME}/getHDF5Tree`,
  async ({ path }, thunkAPI) => {
    try {
      const response = await axios.get(`${BASE_URL}/hdf5/${path}`)
      return response.data
    } catch (e) {
      return thunkAPI.rejectWithValue(e)
    }
  },
)
