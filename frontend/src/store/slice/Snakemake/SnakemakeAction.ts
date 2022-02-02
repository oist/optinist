import { createAsyncThunk } from '@reduxjs/toolkit'
import { Snakemake_SLICE_NAME, SnakemakeListDTO } from './SnakemakeType'
import axios from 'axios'
import { BASE_URL } from 'const/API'

export const getSnakemakeParams = createAsyncThunk<SnakemakeListDTO, void>(
  `${Snakemake_SLICE_NAME}/getSnakemakeParams`,
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
