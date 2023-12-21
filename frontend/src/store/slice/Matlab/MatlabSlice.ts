import { createSlice } from "@reduxjs/toolkit"

import { getMatlabTree } from "store/slice/Matlab/MatlabAction"
import { MATLAB_SLICE_NAME, MatlabTree } from "store/slice/Matlab/MatlabType"
import { convertToTreeNodeType } from "store/slice/Matlab/MatlabUtils"

const initialState: MatlabTree = {
  isLoading: false,
  tree: [],
}
export const matlabSlice = createSlice({
  name: MATLAB_SLICE_NAME,
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

export default matlabSlice.reducer
