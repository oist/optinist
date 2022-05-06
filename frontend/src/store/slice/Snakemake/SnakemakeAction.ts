import { createAsyncThunk } from '@reduxjs/toolkit'

import { getSnakemakeParamsApi } from 'api/snakemake/Snakemake'
import { ParamDTO } from 'utils/param/ParamType'
import { SNAKEMAKE_SLICE_NAME } from './SnakemakeType'

export const getSnakemakeParams = createAsyncThunk<ParamDTO, void>(
  `${SNAKEMAKE_SLICE_NAME}/getSnakemakeParams`,
  async (_, thunkAPI) => {
    const { rejectWithValue } = thunkAPI
    try {
      const response = getSnakemakeParamsApi()
      return response
    } catch (e) {
      return rejectWithValue(e)
    }
  },
)
