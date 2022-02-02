import { createSlice, PayloadAction } from '@reduxjs/toolkit'
import { getSnakemakeParams } from './SnakemakeAction'
import { SnakemakeType, Snakemake_SLICE_NAME } from './SnakemakeType'
import { convertToSnakemakeListType, getSnakemakeChild } from './SnakemakeUtils'

const initialState: SnakemakeType = {
  SnakemakeList: {},
}

export const SnakemakeSlice = createSlice({
  name: Snakemake_SLICE_NAME,
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
      const targetNode = getSnakemakeChild(state.SnakemakeList, paramPath)
      if (targetNode != null) {
        targetNode.value = newValue
      }
    },
  },
  extraReducers: (builder) => {
    builder.addCase(getSnakemakeParams.fulfilled, (state, action) => {
      state.SnakemakeList = convertToSnakemakeListType(action.payload)
    })
  },
})

export const { updateParam } = SnakemakeSlice.actions

export default SnakemakeSlice.reducer
