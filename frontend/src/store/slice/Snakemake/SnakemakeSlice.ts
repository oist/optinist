import { createSlice, PayloadAction } from "@reduxjs/toolkit"

import { clearFlowElements } from "store/slice/FlowElement/FlowElementSlice"
import { getSnakemakeParams } from "store/slice/Snakemake/SnakemakeAction"
import {
  SnakemakeType,
  SNAKEMAKE_SLICE_NAME,
} from "store/slice/Snakemake/SnakemakeType"
import { convertToParamMap, getChildParam } from "utils/param/ParamUtils"

const initialState: SnakemakeType = {
  params: {},
}

export const SnakemakeSlice = createSlice({
  name: SNAKEMAKE_SLICE_NAME,
  initialState,
  reducers: {
    updateParam: (
      state,
      action: PayloadAction<{
        path: string
        newValue: unknown
      }>,
    ) => {
      const { path, newValue } = action.payload
      const target = getChildParam(path, state.params)
      if (target != null) {
        target.value = newValue
      }
    },
  },
  extraReducers: (builder) => {
    builder
      .addCase(getSnakemakeParams.fulfilled, (state, action) => {
        state.params = convertToParamMap(action.payload)
      })
      .addCase(clearFlowElements, () => initialState)
  },
})

export const { updateParam } = SnakemakeSlice.actions

export default SnakemakeSlice.reducer
