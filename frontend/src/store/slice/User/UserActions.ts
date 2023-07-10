import { createAsyncThunk } from '@reduxjs/toolkit'
import { USER_SLICE_NAME } from './UserType'
import { deleteMeApi, getMeApi, updateMeApi } from 'api/users/UsersMe'
import { UpdateUserDTO } from 'api/users/UsersApiDTO'
import { LoginDTO, loginApi } from 'api/auth/Auth'

export const login = createAsyncThunk(
  `${USER_SLICE_NAME}/login`,
  async (data: LoginDTO, thunkAPI) => {
    try {
      const responseData = await loginApi(data)
      return responseData
    } catch (e) {
      return thunkAPI.rejectWithValue(e)
    }
  },
)

export const getMe = createAsyncThunk(
  `${USER_SLICE_NAME}/getMe`,
  async (_, thunkAPI) => {
    try {
      const responseData = await getMeApi()
      return responseData
    } catch (e) {
      return thunkAPI.rejectWithValue(e)
    }
  },
)

export const updateMe = createAsyncThunk(
  `${USER_SLICE_NAME}/updateMe`,
  async (data: UpdateUserDTO, thunkAPI) => {
    try {
      const responseData = await updateMeApi(data)
      return responseData
    } catch (e) {
      return thunkAPI.rejectWithValue(e)
    }
  },
)

export const deleteMe = createAsyncThunk(
  `${USER_SLICE_NAME}/deleteMe`,
  async (_, thunkAPI) => {
    try {
      const responseData = await deleteMeApi()
      return responseData
    } catch (e) {
      return thunkAPI.rejectWithValue(e)
    }
  },
)
