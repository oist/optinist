import { createAsyncThunk } from '@reduxjs/toolkit'
import { NWB_SLICE_NAME } from './NWBType'
import axios from 'axios'
import { BASE_URL } from 'const/API'
import { ParamDTO } from 'store/utils/param/ParamType'

export const getNWBParams = createAsyncThunk<ParamDTO, void>(
  `${NWB_SLICE_NAME}/getNWBParams`,
  async (_, thunkAPI) => {
    const { rejectWithValue } = thunkAPI
    try {
      const response = await axios.get(`${BASE_URL}/nwb`)
      return response.data
    } catch (e) {
      return rejectWithValue(e)
    }
  },
)
