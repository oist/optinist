import { createAsyncThunk } from '@reduxjs/toolkit'
import { getHDF5TreeApi, HDF5TreeDTO } from 'api/hdf5/HDF5'

import { HDF5_SLICE_NAME } from './HDF5Type'

export const getHDF5Tree = createAsyncThunk<
  HDF5TreeDTO[],
  { path: string; workspaceId: number }
>(`${HDF5_SLICE_NAME}/getHDF5Tree`, async ({ path, workspaceId }, thunkAPI) => {
  try {
    const response = await getHDF5TreeApi(path, workspaceId)
    return response
  } catch (e) {
    return thunkAPI.rejectWithValue(e)
  }
})
