import { createAsyncThunk } from "@reduxjs/toolkit"

import { getMatlabTreeApi, MatlabTreeDTO } from "api/matlab/Matlab"
import { HDF5_SLICE_NAME } from "store/slice/HDF5/HDF5Type"

export const getMatlabTree = createAsyncThunk<
  MatlabTreeDTO[],
  { path: string; workspaceId: number }
>(`${HDF5_SLICE_NAME}/getHDF5Tree`, async ({ path, workspaceId }, thunkAPI) => {
  try {
    const response = await getMatlabTreeApi(path, workspaceId)
    return response
  } catch (e) {
    return thunkAPI.rejectWithValue(e)
  }
})
