import { createAsyncThunk } from '@reduxjs/toolkit'
import { ALGORITHM_SLICE_NAME, OutputData } from './AlgorithmType'
import axios from 'axios'
import { AlgoParam } from 'const/NodeData'
import { AlgoOutputDataDTO } from './AlgorithmUtils'

export const getAlgoParams = createAsyncThunk<
  AlgoParam,
  { id: string; algoName: string }
>(`${ALGORITHM_SLICE_NAME}/getAlgoParams`, async ({ algoName }, thunkAPI) => {
  const { rejectWithValue } = thunkAPI
  try {
    const response = await axios.get(
      `http://localhost:8000/api/params/${algoName}`,
    )
    return response.data
  } catch (e) {
    return rejectWithValue(e)
  }
})

export const getAlgoOutputData = createAsyncThunk<
  AlgoOutputDataDTO,
  { id: string; name: string }
>(`${ALGORITHM_SLICE_NAME}/getAlgoOutputData`, async ({ name }, thunkAPI) => {
  try {
    const response = await axios.get(`http://localhost:8000/output/${name}`)
    return response.data
  } catch (e) {
    return thunkAPI.rejectWithValue(e)
  }
})
