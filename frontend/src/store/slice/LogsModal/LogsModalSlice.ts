import { createSlice } from "@reduxjs/toolkit"

type LogsModal = {
  open: boolean
}

const initialState: LogsModal = {
  open: false,
}

export const logsModalSlice = createSlice({
  name: "logsModal",
  initialState,
  reducers: {
    openLogsModal: (state) => {
      state.open = true
    },
    closeLogsModal: (state) => {
      state.open = false
    },
    toggleLogsModal: (state) => {
      state.open = !state.open
    },
  },
})

export const { openLogsModal, closeLogsModal, toggleLogsModal } =
  logsModalSlice.actions

export default logsModalSlice.reducer
