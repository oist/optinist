import { createSlice, PayloadAction } from '@reduxjs/toolkit'
import { getNWBParams } from './NWBAction'
import { NWBType, NWB_SLICE_NAME } from './NWBType'
import { convertToNWBListType } from './NWBUtils'

const initialState: NWBType = {
  params: {},
  nwbList: {},
}

export const nwbSlice = createSlice({
  name: NWB_SLICE_NAME,
  initialState,
  reducers: {
    updateParam: (
      state,
      action: PayloadAction<{
        paramName: string
        newValue: unknown
      }>,
    ) => {
      const { paramName, newValue } = action.payload
      state.params[paramName] = newValue
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
