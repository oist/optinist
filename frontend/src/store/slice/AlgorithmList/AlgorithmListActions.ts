import { createAsyncThunk } from "@reduxjs/toolkit"

import { AlgoListDTO, getAlgoListApi } from "api/algolist/AlgoList"
import { ALGORITHM_LIST_SLICE_NAME } from "store/slice/AlgorithmList/AlgorithmListType"

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
