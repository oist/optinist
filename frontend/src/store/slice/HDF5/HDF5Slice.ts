import { createSlice } from "@reduxjs/toolkit"

import { getHDF5Tree } from "store/slice/HDF5/HDF5Action"
import { HDF5Tree, HDF5_SLICE_NAME } from "store/slice/HDF5/HDF5Type"
import { convertToTreeNodeType } from "store/slice/HDF5/HDF5Utils"

const initialState: HDF5Tree = {
  isLoading: false,
  tree: [],
}
export const HDF5Slice = createSlice({
  name: HDF5_SLICE_NAME,
  initialState,
  reducers: {},
  extraReducers: (builder) => {
    builder
      .addCase(getHDF5Tree.pending, (state) => {
        state.tree = []
        state.isLoading = true
      })
      .addCase(getHDF5Tree.fulfilled, (state, action) => {
        state.tree = convertToTreeNodeType(action.payload)
        state.isLoading = false
      })
  },
})

export default HDF5Slice.reducer
