import { createAsyncThunk } from '@reduxjs/toolkit'

import { getNWBParamsApi } from 'api/nwb/NWB'
import { ParamDTO } from 'utils/param/ParamType'
import { NWB_SLICE_NAME } from './NWBType'

export const getNWBParams = createAsyncThunk<ParamDTO, void>(
  `${NWB_SLICE_NAME}/getNWBParams`,
  async (_, thunkAPI) => {
    const { rejectWithValue } = thunkAPI
    try {
      const response = await getNWBParamsApi()
      return response
    } catch (e) {
      return rejectWithValue(e)
    }
  },
)
