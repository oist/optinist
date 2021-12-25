import { createSlice, PayloadAction } from '@reduxjs/toolkit'
import { getNWBParams } from './NWBAction'
import { NWBType, NWB_SLICE_NAME } from './NWBType'
import { convertToNWBListType, getNWBChild } from './NWBUtils'

const initialState: NWBType = {
  nwbList: {},
}

export const nwbSlice = createSlice({
  name: NWB_SLICE_NAME,
  initialState,
  reducers: {
    updateParam: (
      state,
      action: PayloadAction<{
        paramPath: string
        newValue: unknown
      }>,
    ) => {
      const { paramPath, newValue } = action.payload
      const targetNode = getNWBChild(state.nwbList, paramPath)
      if (targetNode != null) {
        targetNode.value = newValue
      }
    },
  },
  extraReducers: (builder) => {
    builder.addCase(getNWBParams.fulfilled, (state, action) => {
      state.nwbList = convertToNWBListType(action.payload)
    })
  },
})

export const { updateParam } = nwbSlice.actions

export default nwbSlice.reducer
