import { createAsyncThunk } from '@reduxjs/toolkit'
import { ALGORITHM_LIST_SLICE_NAME } from './AlgorithmListType'
import { AlgoListDTO, getAlgoListApi } from 'api/algolist/AlgoList'

export const getAlgoList = createAsyncThunk<AlgoListDTO, void>(
  `${ALGORITHM_LIST_SLICE_NAME}/getAlgoList`,
  async (_, thunkAPI) => {
    const { rejectWithValue } = thunkAPI
    try {
      const response = await getAlgoListApi()
      return response
    } catch (e) {
      return rejectWithValue(e)
    }
  },
)
