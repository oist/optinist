import { createSlice } from "@reduxjs/toolkit"

import { getAlgoList } from "store/slice/AlgorithmList/AlgorithmListActions"
import {
  ALGORITHM_LIST_SLICE_NAME,
  AlgorithmListType,
} from "store/slice/AlgorithmList/AlgorithmListType"
import { convertToAlgoListType } from "store/slice/AlgorithmList/AlgorithmListUtils"
import { getWorkspace } from "store/slice/Workspace/WorkspaceActions"

export const initialState: AlgorithmListType = {
  isLatest: false,
  tree: {},
}

export const algorithmListSlice = createSlice({
  name: ALGORITHM_LIST_SLICE_NAME,
  initialState,
  reducers: {
    resetAlgoListCache: (state) => {
      state.isLatest = false
    },
  },
  extraReducers: (builder) => {
    builder
      .addCase(getAlgoList.fulfilled, (state, action) => {
        state.tree = convertToAlgoListType(action.payload)
        state.isLatest = true
      })
      .addCase(getWorkspace.fulfilled, (state) => {
        // Reset AlgoList cache when workspace is switched
        state.isLatest = false
      })
  },
})

export const { resetAlgoListCache } = algorithmListSlice.actions
export default algorithmListSlice.reducer
