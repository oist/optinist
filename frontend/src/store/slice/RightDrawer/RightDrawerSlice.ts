import { createSlice, PayloadAction } from '@reduxjs/toolkit'

type RightDrawer = {
  open: boolean
  mode: RIGHT_DRAWER_MODE_TYPE
  currendNodeId: string | null
}

export const RIGHT_DRAWER_MODE = {
  NONE: 'none',
  NWB: 'nwb',
  PARAM_FORM: 'paramForm',
} as const

export type RIGHT_DRAWER_MODE_TYPE =
  typeof RIGHT_DRAWER_MODE[keyof typeof RIGHT_DRAWER_MODE]

const initialState: RightDrawer = {
  open: false,
  mode: RIGHT_DRAWER_MODE.NONE,
  currendNodeId: null,
}
export const rightDrawerSlice = createSlice({
  name: 'rightDrawer',
  initialState,
  reducers: {
    openRightDrawer: (state, action: PayloadAction<RIGHT_DRAWER_MODE_TYPE>) => {
      state.open = true
      state.mode = action.payload
    },
    closeRightDrawer: (state) => {
      state.open = false
      state.mode = RIGHT_DRAWER_MODE.NONE
    },
    toggleParamForm: (state, action: PayloadAction<string>) => {
      if (
        state.open &&
        state.mode === RIGHT_DRAWER_MODE.PARAM_FORM &&
        state.currendNodeId === action.payload
      ) {
        state.open = false
        state.mode = RIGHT_DRAWER_MODE.NONE
        state.currendNodeId = null
      } else {
        state.open = true
        state.mode = RIGHT_DRAWER_MODE.PARAM_FORM
        state.currendNodeId = action.payload
      }
    },
    toggleNwb: (state) => {
      if (state.open && state.mode === RIGHT_DRAWER_MODE.NWB) {
        state.open = false
        state.mode = RIGHT_DRAWER_MODE.NONE
      } else {
        state.open = true
        state.mode = RIGHT_DRAWER_MODE.NWB
      }
      state.currendNodeId = null
    },
  },
})

export const { toggleParamForm, toggleNwb, openRightDrawer, closeRightDrawer } =
  rightDrawerSlice.actions

export default rightDrawerSlice.reducer