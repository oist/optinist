import { RootState } from "store/store"

export const selectLogsModalIsOpen = (state: RootState) => state.logsModal.open
