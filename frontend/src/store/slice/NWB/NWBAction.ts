import { createAsyncThunk } from '@reduxjs/toolkit'
import { NWB_SLICE_NAME, NWBListDTO } from './NWBType'
import axios from 'axios'
import { BASE_URL } from 'const/API'

export const getNWBParams = createAsyncThunk<NWBListDTO, void>(
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
