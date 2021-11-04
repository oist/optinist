import { createAsyncThunk } from '@reduxjs/toolkit'
import { ALGORITHM_SLICE_NAME } from './AlgorithmType'
import axios from 'axios'
import { AlgoParam } from 'const/NodeData'
import { BASE_URL } from 'const/API'

export const getAlgoParams = createAsyncThunk<
  AlgoParam,
  { id: string; algoName: string }
>(`${ALGORITHM_SLICE_NAME}/getAlgoParams`, async ({ algoName }, thunkAPI) => {
  const { rejectWithValue } = thunkAPI
  try {
    const response = await axios.get(`${BASE_URL}/api/params/${algoName}`)
    return response.data
  } catch (e) {
    return rejectWithValue(e)
  }
})
