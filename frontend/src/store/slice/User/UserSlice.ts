import { createSlice, isAnyOf } from '@reduxjs/toolkit'
import { USER_SLICE_NAME } from './UserType'
import { User } from './UserType'
import { deleteMe, getMe, login, updateMe } from './UserActions'
import {
  removeExToken,
  removeToken,
  saveExToken,
  saveRefreshToken,
  saveToken,
} from 'utils/auth/AuthUtils'

const initialState: User = { currentUser: undefined }

export const userSlice = createSlice({
  name: USER_SLICE_NAME,
  initialState,
  reducers: {},
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
      .addCase(updateMe.fulfilled, (state, action) => {
        state.currentUser = action.payload
      })
      .addMatcher(
        isAnyOf(login.rejected, getMe.rejected, deleteMe.fulfilled),
        (state, _) => {
          removeToken()
          removeExToken()
          state = initialState
        },
      )
  },
})

export default userSlice.reducer
