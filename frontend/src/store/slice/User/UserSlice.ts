import { createSlice, isAnyOf } from '@reduxjs/toolkit'
import { USER_SLICE_NAME } from './UserType'
import { User } from './UserType'
import {
  deleteMe,
  getListUser,
  getListSearch,
  getMe,
  login,
  updateMe,
  updateMePassword,
  deleteUser,
  createUser,
  updateUser,
} from './UserActions'
import {
  removeExToken,
  removeToken,
  saveExToken,
  saveRefreshToken,
  saveToken,
} from 'utils/auth/AuthUtils'

const initialState: User = {
  currentUser: undefined,
  listUserSearch: undefined,
  listUser: undefined,
  loading: false,
}

export const userSlice = createSlice({
  name: USER_SLICE_NAME,
  initialState,
  reducers: {
    logout: (state) => {
      removeToken()
      removeExToken()
      state = initialState
    },
    resetUserSearch: (state) => {
      state.listUserSearch = []
    },
  },
  extraReducers: (builder) => {
    builder
      .addCase(login.fulfilled, (_, action) => {
        saveToken(action.payload.access_token)
        saveRefreshToken(action.payload.refresh_token)
        saveExToken(action.payload.ex_token)
      })
      .addCase(getMe.fulfilled, (state, action) => {
        state.currentUser = action.payload
      })
      .addCase(getListUser.fulfilled, (state, action) => {
        state.listUser = action.payload
        state.loading = false
      })
      .addCase(getListSearch.fulfilled, (state, action) => {
        state.listUserSearch = action.payload
      })
      .addMatcher(
        isAnyOf(
          getListSearch.rejected,
          createUser.rejected,
          getListUser.rejected,
          updateUser.rejected,
          updateMePassword.rejected,
          updateMePassword.fulfilled,
          deleteUser.fulfilled,
          deleteUser.rejected,
          deleteMe.rejected,
          deleteMe.fulfilled,
          updateMe.rejected,
          updateMe.fulfilled,
          createUser.fulfilled,
        ),
        (state) => {
          state.loading = false
        },
      )
      .addMatcher(
        isAnyOf(
          getListUser.pending,
          deleteUser.pending,
          createUser.pending,
          updateMe.pending,
          deleteMe.pending,
          updateUser.pending,
          updateMePassword.pending,
        ),
        (state) => {
          state.loading = true
        },
      )
      .addMatcher(
        isAnyOf(login.rejected, getMe.rejected, deleteMe.fulfilled),
        (state) => {
          removeToken()
          removeExToken()
          state = initialState
        },
      )
  },
})

export const { logout, resetUserSearch } = userSlice.actions
export default userSlice.reducer
