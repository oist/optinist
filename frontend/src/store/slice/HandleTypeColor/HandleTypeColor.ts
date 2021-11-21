import { createSlice, PayloadAction } from '@reduxjs/toolkit'

export type HandleTypeColor = {
  colorMap: { [type: string]: string }
  /**
   * HANDLE_COLOR_PRESET_MAPのkey
   */
  nextKey: number
}

const initialState: HandleTypeColor = { colorMap: {}, nextKey: 0 }

const SLICE_NAME = 'handleTypeColor'

export const handleTypeColorSlice = createSlice({
  name: SLICE_NAME,
  initialState,
  reducers: {
    addColor: (
      state,
      action: PayloadAction<{ type: string; color: string }>,
    ) => {
      state.colorMap[action.payload.type] = action.payload.color
      state.nextKey++
    },
  },
})

export const { addColor } = handleTypeColorSlice.actions

export default handleTypeColorSlice.reducer
