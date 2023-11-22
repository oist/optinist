import { createSlice } from "@reduxjs/toolkit"

import { getModeStandalone } from "store/slice/Standalone/StandaloneActions"
import { ModeType } from "store/slice/Standalone/StandaloneType"

const initialState: ModeType = {
  mode: undefined,
  loading: false,
}

export const modeStandaloneSlice = createSlice({
  name: "modeStandalone",
  initialState,
  reducers: {},
  extraReducers(build) {
    build
      .addCase(getModeStandalone.pending, (state) => {
        state.loading = true
      })
      .addCase(getModeStandalone.rejected, (state) => {
        state.loading = false
      })
      .addCase(getModeStandalone.fulfilled, (state, action) => {
        state.mode = action.payload
        state.loading = false
      })
  },
})

export default modeStandaloneSlice.reducer
