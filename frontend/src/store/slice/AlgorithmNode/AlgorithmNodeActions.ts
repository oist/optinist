import { createAsyncThunk } from '@reduxjs/toolkit'
import { getAlgoParamsApi } from 'api/params/Params'
import { ParamDTO } from 'utils/param/ParamType'
import { ALGORITHM_NODE_SLICE_NAME } from './AlgorithmNodeType'

export const getAlgoParams = createAsyncThunk<
  ParamDTO,
  { nodeId: string; algoName: string }
>(
  `${ALGORITHM_NODE_SLICE_NAME}/getAlgoParams`,
  async ({ algoName }, thunkAPI) => {
    const { rejectWithValue } = thunkAPI
    try {
      const response = await getAlgoParamsApi(algoName)
      return response
    } catch (e) {
      return rejectWithValue(e)
    }
  },
)
