import { createSlice, PayloadAction } from '@reduxjs/toolkit'
import { HANDLE_COLOR_PRESET_MAP } from 'const/HandleColor'

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
    addColor: (state, action: PayloadAction<string>) => {
      const nextColor =
        HANDLE_COLOR_PRESET_MAP.get(state.nextKey) ??
        '#' + Math.floor(Math.random() * 0xffffff).toString(16)
      state.colorMap[action.payload] = nextColor
      state.nextKey++
    },
  },
})

export const { addColor } = handleTypeColorSlice.actions

export default handleTypeColorSlice.reducer
