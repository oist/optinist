import { createAsyncThunk } from '@reduxjs/toolkit'
import axios from 'axios'
import { AlgoListDTO, ALGORITHM_LIST_SLICE_NAME } from './AlgorithmListType'
import { BASE_URL } from 'const/API'

export const getAlgoList = createAsyncThunk<AlgoListDTO, void>(
  `${ALGORITHM_LIST_SLICE_NAME}/getAlgoList`,
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
