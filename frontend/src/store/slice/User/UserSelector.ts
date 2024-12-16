import { ROLE } from "@types"
import { RootState } from "store/store"

export const selectCurrentUser = (state: RootState) => state.user.currentUser
export const selectCurrentUserId = (state: RootState) =>
  selectCurrentUser(state)?.id
export const selectListUser = (state: RootState) => state.user.listUser
export const selectLoading = (state: RootState) => state.user.loading
export const selectCurrentUserUid = (state: RootState) =>
  selectCurrentUser(state)?.uid
export const selectCurrentUserEmail = (state: RootState) =>
  selectCurrentUser(state)?.email
export const selectListUserSearch = (state: RootState) =>
  state.user.listUserSearch
export const isAdmin = (state: RootState) => {
  return state.user && ROLE.ADMIN === state.user.currentUser?.role_id
}
