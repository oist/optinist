import { createSlice } from "@reduxjs/toolkit"

import { HDF5Tree, HDF5_SLICE_NAME } from "store/slice/HDF5/HDF5Type"
import { convertToTreeNodeType } from "store/slice/HDF5/HDF5Utils"
import { getMatlabTree } from "store/slice/Matlab/MatlabAction"

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
      .addCase(getMatlabTree.pending, (state) => {
        state.tree = []
        state.isLoading = true
      })
      .addCase(getMatlabTree.fulfilled, (state, action) => {
        state.tree = convertToTreeNodeType(action.payload)
        state.isLoading = false
      })
  },
})

export default HDF5Slice.reducer
