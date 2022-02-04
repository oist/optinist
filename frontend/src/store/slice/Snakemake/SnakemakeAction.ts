import { createAsyncThunk } from '@reduxjs/toolkit'
import { SNAKEMAKE_SLICE_NAME } from './SnakemakeType'
import axios from 'axios'
import { BASE_URL } from 'const/API'
import { ParamDTO } from 'store/utils/param/ParamType'

export const getSnakemakeParams = createAsyncThunk<ParamDTO, void>(
  `${SNAKEMAKE_SLICE_NAME}/getSnakemakeParams`,
  async (_, thunkAPI) => {
    const { rejectWithValue } = thunkAPI
    try {
      const response = await axios.get(`${BASE_URL}/snakemake`)
      return response.data
    } catch (e) {
      return rejectWithValue(e)
    }
  },
)
