import { createAsyncThunk } from '@reduxjs/toolkit'
import axios from 'axios'
import { BASE_URL } from 'const/API'
import {
  ImageData,
  TimeSeriesData,
  DISPLAY_DATA_SLICE_NAME,
  HeatMapData,
  TableData,
  RoiData,
  ScatterData,
} from './DisplayDataType'

export const getTimeSeriesData = createAsyncThunk<
  { data: TimeSeriesData },
  { path: string; index: number }
>(
  `${DISPLAY_DATA_SLICE_NAME}/getTimeSeriesData`,
  async ({ path, index }, thunkAPI) => {
    try {
      const response = await axios.get(`${BASE_URL}/outputs/timedata/${path}`, {
        params: {
          index: index,
        },
      })
      return response.data
    } catch (e) {
      return thunkAPI.rejectWithValue(e)
    }
  },
)

export const getHeatMapData = createAsyncThunk<
  { data: HeatMapData },
  { path: string }
>(`${DISPLAY_DATA_SLICE_NAME}/getHeatMapData`, async ({ path }, thunkAPI) => {
  try {
    const response = await axios.get(`${BASE_URL}/outputs/data/${path}`)
    return response.data
  } catch (e) {
    return thunkAPI.rejectWithValue(e)
  }
})

export const getImageData = createAsyncThunk<
  { data: ImageData },
  { path: string; startIndex?: number; endIndex?: number }
>(
  `${DISPLAY_DATA_SLICE_NAME}/getImageData`,
  async ({ path, startIndex, endIndex }, thunkAPI) => {
    try {
      const response = await axios.get(`${BASE_URL}/outputs/image/${path}`, {
        params: {
          start_index: startIndex,
          end_index: endIndex,
        },
      })
      return response.data
    } catch (e) {
      return thunkAPI.rejectWithValue(e)
    }
  },
)

export const getTableData = createAsyncThunk<
  {
    columns: string[]
    data: TableData
  },
  { path: string }
>(`${DISPLAY_DATA_SLICE_NAME}/getTableData`, async ({ path }, thunkAPI) => {
  try {
    const response = await axios.get(`${BASE_URL}/outputs/csv/${path}`)
    return response.data.data
  } catch (e) {
    return thunkAPI.rejectWithValue(e)
  }
})

export const getRoiData = createAsyncThunk<{ data: RoiData }, { path: string }>(
  `${DISPLAY_DATA_SLICE_NAME}/getRoiData`,
  async ({ path }, thunkAPI) => {
    try {
      const response = await axios.get(`${BASE_URL}/outputs/image/${path}`, {})
      return response.data
    } catch (e) {
      return thunkAPI.rejectWithValue(e)
    }
  },
)

export const getScatterData = createAsyncThunk<
  { data: ScatterData },
  { path: string }
>(`${DISPLAY_DATA_SLICE_NAME}/getScatterData`, async ({ path }, thunkAPI) => {
  try {
    const response = await axios.get(`${BASE_URL}/outputs/data/${path}`, {})
    return response.data
  } catch (e) {
    return thunkAPI.rejectWithValue(e)
  }
})
