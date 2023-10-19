import { createAsyncThunk } from "@reduxjs/toolkit"
import { USER_SLICE_NAME } from "./UserType"
import {
  deleteMeApi,
  getMeApi,
  updateMeApi,
  updateMePasswordApi,
} from "api/users/UsersMe"
import {
  AddUserDTO,
  ListUsersQueryDTO,
  UpdateUserDTO,
  UpdateUserPasswordDTO,
  UserDTO,
} from "api/users/UsersApiDTO"
import { LoginDTO, loginApi } from "api/auth/Auth"
import {
  deleteUserApi,
  createUserApi,
  listUsersApi,
  updateUserApi,
} from "api/users/UsersAdmin"
import { getListSearchApi } from "api/users/UsersAdmin"

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
      await thunkAPI.dispatch(getMe())
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

export const updateMePassword = createAsyncThunk(
  `${USER_SLICE_NAME}/updateMePassword`,
  async (data: UpdateUserPasswordDTO, thunkAPI) => {
    try {
      return await updateMePasswordApi(data)
    } catch (e) {
      return thunkAPI.rejectWithValue(e)
    }
  },
)

export const getListUser = createAsyncThunk(
  `${USER_SLICE_NAME}/getListUser`,
  async (params: ListUsersQueryDTO, thunkAPI) => {
    try {
      return await listUsersApi(params)
    } catch (e) {
      return thunkAPI.rejectWithValue(e)
    }
  },
)

export const getListSearch = createAsyncThunk<
  UserDTO[],
  { keyword: string | null }
>(`${USER_SLICE_NAME}/getListSearch`, async (params, thunkAPI) => {
  try {
    const responseData = await getListSearchApi(params)
    return responseData
  } catch (e) {
    return thunkAPI.rejectWithValue(e)
  }
})

export const createUser = createAsyncThunk<
  UserDTO,
  { data: AddUserDTO; params: ListUsersQueryDTO }
>(`${USER_SLICE_NAME}/createUser`, async (props, thunkAPI) => {
  try {
    const { dispatch } = thunkAPI
    const responseData = await createUserApi(props.data)
    await dispatch(getListUser(props.params))
    return responseData
  } catch (e) {
    return thunkAPI.rejectWithValue(e)
  }
})

export const updateUser = createAsyncThunk<
  UserDTO,
  { id: number; data: UpdateUserDTO; params: ListUsersQueryDTO }
>(`${USER_SLICE_NAME}/updateUser`, async (props, thunkAPI) => {
  const { dispatch } = thunkAPI
  try {
    const responseData = await updateUserApi(props.id, props.data)
    await dispatch(getListUser(props.params))
    return responseData
  } catch (e) {
    return thunkAPI.rejectWithValue(e)
  }
})

export const deleteUser = createAsyncThunk<
  string,
  { id: number; params: ListUsersQueryDTO }
>(`${USER_SLICE_NAME}/deleteUser`, async (data, thunkAPI) => {
  const { dispatch } = thunkAPI
  try {
    const responseData = await deleteUserApi(data.id)
    await dispatch(getListUser(data.params))
    return responseData
  } catch (e) {
    return thunkAPI.rejectWithValue(e)
  }
})
