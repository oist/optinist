import { createSlice } from '@reduxjs/toolkit'
import {
  ALGORITHM_LIST_SLICE_NAME,
  AlgorithmListType,
} from './AlgorithmListType'
import { getAlgoList } from './AlgorithmListActions'
import { convertToAlgoListType } from './AlgorithmListUtils'

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
