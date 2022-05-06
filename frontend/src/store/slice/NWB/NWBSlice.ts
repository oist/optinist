import { createSlice, PayloadAction } from '@reduxjs/toolkit'
import { getNWBParams } from './NWBAction'
import { NWBType, NWB_SLICE_NAME } from './NWBType'
import { convertToParamMap, getChildParam } from 'utils/param/ParamUtils'

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
    builder.addCase(getNWBParams.fulfilled, (state, action) => {
      state.params = convertToParamMap(action.payload)
    })
  },
})

export const { updateParam } = nwbSlice.actions

export default nwbSlice.reducer
