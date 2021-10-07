import axios from 'axios'
import { OUTPUT_SLICE_NAME, OutputData } from './OutputType'
import { createAsyncThunk } from '@reduxjs/toolkit'

export const getOutputData = createAsyncThunk<OutputData[], { id: string }>(
  `${OUTPUT_SLICE_NAME}/getOutputData`,
  async ({ id }, thunkAPI) => {
    try {
      const response = await axios.get(`http://localhost:8000/output`)
      return response.data.data
    } catch (e) {
      return thunkAPI.rejectWithValue(e)
    }
  },
)
