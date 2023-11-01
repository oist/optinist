import { createSlice, PayloadAction } from "@reduxjs/toolkit"

import { clearFlowElements } from "store/slice/FlowElement/FlowElementSlice"
import { getNWBParams } from "store/slice/NWB/NWBAction"
import { NWBType, NWB_SLICE_NAME } from "store/slice/NWB/NWBType"
import { convertToParamMap, getChildParam } from "utils/param/ParamUtils"
const initialState: NWBType = {
  params: {},
}

export const nwbSlice = createSlice({
  name: NWB_SLICE_NAME,
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
      const targetNode = getChildParam(path, state.params)
      if (targetNode != null) {
        targetNode.value = newValue
      }
    },
  },
  extraReducers: (builder) => {
    builder
      .addCase(getNWBParams.fulfilled, (state, action) => {
        state.params = convertToParamMap(action.payload)
      })
      .addCase(clearFlowElements, () => initialState)
  },
})

export const { updateParam } = nwbSlice.actions

export default nwbSlice.reducer
