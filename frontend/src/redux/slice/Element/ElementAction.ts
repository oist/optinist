import { createAsyncThunk } from '@reduxjs/toolkit'
import axios from 'axios'
import { Param } from './ElementType'

export const getAlgoParams = createAsyncThunk<
  Param,
  { id: string; algoName: string }
>('element/getAlgoParams', async ({ algoName }, thunkAPI) => {
  const { rejectWithValue } = thunkAPI
  try {
    const response = await axios.get(`http://localhost:8000/params/${algoName}`)
    return response.data
  } catch (e) {
    return rejectWithValue(e)
  }
})
