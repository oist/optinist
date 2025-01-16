import { RootState } from "store/store"

export const selectModeStandalone = (state: RootState) => state.mode.mode
export const selectLoading = (state: RootState) => state.mode.loading
