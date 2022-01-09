import { RootState } from 'store/store'

export const selectRightDrawerIsOpen = (state: RootState) =>
  state.rightDrawer.open

export const selectRightDrawerMode = (state: RootState) =>
  state.rightDrawer.mode

export const selectRightDrawerCurrentNodeId = (state: RootState) =>
  state.rightDrawer.currendNodeId
