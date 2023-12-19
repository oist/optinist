import { createSlice } from "@reduxjs/toolkit"

import { getAlgoList } from "store/slice/AlgorithmList/AlgorithmListActions"
import {
  ALGORITHM_LIST_SLICE_NAME,
  AlgorithmListType,
} from "store/slice/AlgorithmList/AlgorithmListType"
import { convertToAlgoListType } from "store/slice/AlgorithmList/AlgorithmListUtils"

export const initialState: AlgorithmListType = {
  isLatest: false,
  tree: {},
}

export const algorithmListSlice = createSlice({
  name: ALGORITHM_LIST_SLICE_NAME,
  initialState,
  reducers: {},
  extraReducers: (builder) => {
    builder.addCase(getAlgoList.fulfilled, (state, action) => {
      state.tree = convertToAlgoListType(action.payload)
      state.isLatest = true
    })
  },
})

export default algorithmListSlice.reducer
