import { createAsyncThunk } from '@reduxjs/toolkit'
import { ALGORITHM_SLICE_NAME, OutputData } from './AlgorithmType'
import axios from 'axios'
import { AlgoParam } from 'const/NodeData'

export const getAlgoParams = createAsyncThunk<
  AlgoParam,
  { id: string; algoName: string }
>(`${ALGORITHM_SLICE_NAME}/getAlgoParams`, async ({ algoName }, thunkAPI) => {
  const { rejectWithValue } = thunkAPI
  try {
    const response = await axios.get(`http://localhost:8000/params/${algoName}`)
    return response.data
  } catch (e) {
    return rejectWithValue(e)
  }
})

export const getAlgoOutputData = createAsyncThunk<OutputData[], { id: string }>(
  `${ALGORITHM_SLICE_NAME}/getAlgoOutputData`,
  async ({ id }, thunkAPI) => {
    try {
      const response = await axios.get(`http://localhost:8000/output`)
      return response.data.data
    } catch (e) {
      return thunkAPI.rejectWithValue(e)
    }
  },
)
