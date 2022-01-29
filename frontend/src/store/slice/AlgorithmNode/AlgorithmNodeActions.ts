import { createAsyncThunk } from '@reduxjs/toolkit'
import { ALGORITHM_NODE_SLICE_NAME, ParamDTO } from './AlgorithmNodeType'
import axios from 'axios'
import { BASE_URL } from 'const/API'

export const getAlgoParams = createAsyncThunk<
  ParamDTO,
  { nodeId: string; algoName: string }
>(
  `${ALGORITHM_NODE_SLICE_NAME}/getAlgoParams`,
  async ({ algoName }, thunkAPI) => {
    const { rejectWithValue } = thunkAPI
    try {
      const response = await axios.get(`${BASE_URL}/params/${algoName}`)
      return response.data
    } catch (e) {
      return rejectWithValue(e)
    }
  },
)
