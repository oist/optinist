import { RootState } from "store/store"

export const selectModeStandalone = (state: RootState) => state.mode.mode
export const selectLoadingMode = (state: RootState) => state.mode.loading
