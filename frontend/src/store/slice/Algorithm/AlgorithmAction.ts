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

type AlgoListDTO = {
  [name: string]: {
    args: string[]
    returns: string[]
  } // | { children: AlgoListDTO }
}

export const getAlgoList = createAsyncThunk<AlgoListDTO, void>(
  `${ALGORITHM_SLICE_NAME}/getAlgoList`,
  async (_, thunkAPI) => {
    const { rejectWithValue } = thunkAPI
    try {
      const response = await axios.get(`${BASE_URL}/api/algolist`)
      return response.data
    } catch (e) {
      return rejectWithValue(e)
    }
  },
)
