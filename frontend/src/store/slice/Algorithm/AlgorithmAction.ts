import { createAsyncThunk, createAction } from '@reduxjs/toolkit'
import { AlgoListDTO, ALGORITHM_SLICE_NAME } from './AlgorithmType'
import axios from 'axios'
import { AlgoParam } from 'const/NodeData'
import { BASE_URL } from 'const/API'
import { OutputPathsDTO } from 'api/Run/Run'

export const getAlgoParams = createAsyncThunk<
  AlgoParam,
  { id: string; algoName: string }
>(`${ALGORITHM_SLICE_NAME}/getAlgoParams`, async ({ algoName }, thunkAPI) => {
  const { rejectWithValue } = thunkAPI
  try {
    const response = await axios.get(`${BASE_URL}/params/${algoName}`)
    return response.data
  } catch (e) {
    return rejectWithValue(e)
  }
})

export const getAlgoList = createAsyncThunk<AlgoListDTO, void>(
  `${ALGORITHM_SLICE_NAME}/getAlgoList`,
  async (_, thunkAPI) => {
    const { rejectWithValue } = thunkAPI
    try {
      const response = await axios.get(`${BASE_URL}/algolist`)
      return response.data
    } catch (e) {
      return rejectWithValue(e)
    }
  },
)

export const reflectRunPipelineResult = createAction<{
  dto: OutputPathsDTO
  error?: { name: string; message: string } | null
}>(`${ALGORITHM_SLICE_NAME}/reflectRunPipelineResult`)
