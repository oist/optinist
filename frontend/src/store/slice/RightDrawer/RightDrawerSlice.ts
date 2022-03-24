import { createSlice, PayloadAction } from '@reduxjs/toolkit'
import { importExperimentByUid } from '../Experiments/ExperimentsActions'
import {
  deleteFlowElements,
  deleteFlowElementsById,
} from '../FlowElement/FlowElementSlice'

type RightDrawer = {
  open: boolean
  mode: RIGHT_DRAWER_MODE_TYPE
  currendNodeId: string | null
}

export const RIGHT_DRAWER_MODE = {
  NONE: 'none',
  NWB: 'nwb',
  PARAM_FORM: 'paramForm',
  SNAKEMAKE: 'snakemake',
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
    toggleSnakemake: (state) => {
      if (state.open && state.mode === RIGHT_DRAWER_MODE.SNAKEMAKE) {
        state.open = false
        state.mode = RIGHT_DRAWER_MODE.NONE
      } else {
        state.open = true
        state.mode = RIGHT_DRAWER_MODE.SNAKEMAKE
      }
      state.currendNodeId = null
    },
  },
  extraReducers: (builder) => {
    builder
      .addCase(deleteFlowElements, (state, action) => {
        if (
          action.payload.findIndex((elm) => elm.id === state.currendNodeId) > 0
        ) {
          state.currendNodeId = null
        }
      })
      .addCase(deleteFlowElementsById, (state, action) => {
        if (action.payload === state.currendNodeId) {
          state.currendNodeId = null
        }
      })
      .addCase(importExperimentByUid.fulfilled, () => {
        return initialState
      })
  },
})

export const {
  toggleParamForm,
  toggleNwb,
  toggleSnakemake,
  openRightDrawer,
  closeRightDrawer,
} = rightDrawerSlice.actions

export default rightDrawerSlice.reducer
