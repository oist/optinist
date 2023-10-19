import { createAsyncThunk } from "@reduxjs/toolkit"

import { getAlgoParamsApi } from "api/params/Params"
import { ALGORITHM_NODE_SLICE_NAME } from "store/slice/AlgorithmNode/AlgorithmNodeType"
import { ParamDTO } from "utils/param/ParamType"


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
