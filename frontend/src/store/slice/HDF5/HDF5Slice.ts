import { createSlice } from '@reduxjs/toolkit'

import { getHDF5Tree } from './HDF5Action'
import { HDF5Tree, HDF5_SLICE_NAME } from './HDF5Type'
import { convertToTreeNodeType } from './HDF5Utils'

const initialState: HDF5Tree = {
  isLatest: false,
  isLoading: false,
  tree: [],
}
export const HDF5Slice = createSlice({
  name: HDF5_SLICE_NAME,
  initialState,
  reducers: {},
  extraReducers: (builder) => {
    builder
      .addCase(getHDF5Tree.pending, (state, action) => {
        state = {
          isLoading: true,
          isLatest: false,
          tree: [],
        }
      })
      .addCase(getHDF5Tree.fulfilled, (state, action) => {
        state.tree = convertToTreeNodeType(action.payload)
        state.isLatest = true
        state.isLoading = false
      })
  },
})

export default HDF5Slice.reducer
