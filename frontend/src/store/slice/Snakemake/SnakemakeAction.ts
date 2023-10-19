import { createAsyncThunk } from "@reduxjs/toolkit"

import { getSnakemakeParamsApi } from "api/snakemake/Snakemake"
import { SNAKEMAKE_SLICE_NAME } from "store/slice/Snakemake/SnakemakeType"
import { ParamDTO } from "utils/param/ParamType"


export const getSnakemakeParams = createAsyncThunk<ParamDTO, void>(
  `${SNAKEMAKE_SLICE_NAME}/getSnakemakeParams`,
  async (_, thunkAPI) => {
    const { rejectWithValue } = thunkAPI
    try {
      const response = getSnakemakeParamsApi()
      return response
    } catch (e) {
      return rejectWithValue(e)
    }
  },
)
