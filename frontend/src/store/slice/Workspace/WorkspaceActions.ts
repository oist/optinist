import { createAsyncThunk } from "@reduxjs/toolkit"
import {
  WorkspacePostDataDTO,
  delWorkspaceApi,
  getWorkspaceApi,
  getWorkspacesApi,
  postWorkspaceApi,
  putWorkspaceApi,
  getListUserShareWorkspaceApi,
  postListUserShareWorkspaceApi,
} from "api/workspace"
import {
  ItemsWorkspace,
  WorkspaceDataDTO,
  WorkspaceParams,
  ListUserShareWorkspaceDTO,
  WORKSPACE_SLICE_NAME,
} from "./WorkspaceType"

export const getWorkspace = createAsyncThunk<ItemsWorkspace, { id: number }>(
  `${WORKSPACE_SLICE_NAME}/getWorkspace`,
  async (data, thunkAPI) => {
    const { rejectWithValue } = thunkAPI
    try {
      const response = await getWorkspaceApi(data.id)
      return response
    } catch (e) {
      return rejectWithValue(e)
    }
  },
)

export const getWorkspaceList = createAsyncThunk<
  WorkspaceDataDTO,
  { [key: string]: number }
>(`${WORKSPACE_SLICE_NAME}/getWorkspaceList`, async (params, thunkAPI) => {
  const { rejectWithValue } = thunkAPI
  try {
    const response = await getWorkspacesApi(params)
    return response
  } catch (e) {
    return rejectWithValue(e)
  }
})

export const delWorkspace = createAsyncThunk<boolean, WorkspaceParams>(
  `${WORKSPACE_SLICE_NAME}/delWorkspaceList`,
  async (data, thunkAPI) => {
    const { rejectWithValue, dispatch } = thunkAPI
    try {
      const response = await delWorkspaceApi(Number(data.id))
      await dispatch(getWorkspaceList(data.params as { [key: string]: number }))
      return response
    } catch (e) {
      return rejectWithValue(e)
    }
  },
)

export const postWorkspace = createAsyncThunk<
  ItemsWorkspace,
  WorkspacePostDataDTO
>(`${WORKSPACE_SLICE_NAME}/postWorkspaceList`, async (data, thunkAPI) => {
  const { rejectWithValue } = thunkAPI
  try {
    const response = await postWorkspaceApi(data)
    return response
  } catch (e) {
    return rejectWithValue(e)
  }
})

export const putWorkspace = createAsyncThunk<
  ItemsWorkspace,
  WorkspacePostDataDTO
>(`${WORKSPACE_SLICE_NAME}/putWorkspaceList`, async (data, thunkAPI) => {
  const { rejectWithValue } = thunkAPI
  try {
    const response = await putWorkspaceApi(data)
    return response
  } catch (e) {
    return rejectWithValue(e)
  }
})

export const getListUserShareWorkSpaces = createAsyncThunk<
  ListUserShareWorkspaceDTO,
  { id: number }
>(
  `${WORKSPACE_SLICE_NAME}/getListUserShareWorkSpaces`,
  async (params, thunkAPI) => {
    const { rejectWithValue } = thunkAPI
    try {
      const response = await getListUserShareWorkspaceApi(params.id)
      return response
    } catch (e) {
      return rejectWithValue(e)
    }
  },
)

export const postListUserShareWorkspaces = createAsyncThunk<
  boolean,
  {
    id: number
    data: { user_ids: number[] }
  }
>(
  `${WORKSPACE_SLICE_NAME}/postListUserShareWorkspaces`,
  async (params, thunkAPI) => {
    const { rejectWithValue } = thunkAPI
    try {
      const response = await postListUserShareWorkspaceApi(
        params.id,
        params.data,
      )
      return response
    } catch (e) {
      return rejectWithValue(e)
    }
  },
)
