import { createSlice, PayloadAction } from '@reduxjs/toolkit'
import { convertToParamMap, getChildParam } from 'utils/param/ParamUtils'
import { getSnakemakeParams } from './SnakemakeAction'
import { SnakemakeType, SNAKEMAKE_SLICE_NAME } from './SnakemakeType'

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
    builder.addCase(getSnakemakeParams.fulfilled, (state, action) => {
      state.params = convertToParamMap(action.payload)
    })
  },
})

export const { updateParam } = SnakemakeSlice.actions

export default SnakemakeSlice.reducer
