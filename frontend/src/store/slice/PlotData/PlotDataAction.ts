import { createAsyncThunk } from '@reduxjs/toolkit'
import axios from 'axios'

import { BASE_URL } from 'const/API'
import {
  ImageData,
  TimeSeriesData,
  PLOT_DATA_SLICE_NAME,
  HeatMapData,
} from './PlotDataType'

export const getTimeSeriesData = createAsyncThunk<
  { data: TimeSeriesData },
  { path: string }
>(`${PLOT_DATA_SLICE_NAME}/getTimeSeriesData`, async ({ path }, thunkAPI) => {
  try {
    const response = await axios.get(`${BASE_URL}/outputs/data/${path}`)
    return response.data
  } catch (e) {
    return thunkAPI.rejectWithValue(e)
  }
})

export const getHeatMapData = createAsyncThunk<
  { data: HeatMapData },
  { path: string }
>(`${PLOT_DATA_SLICE_NAME}/getHeatMapData`, async ({ path }, thunkAPI) => {
  try {
    const response = await axios.get(`${BASE_URL}/outputs/data/${path}`)
    return response.data
  } catch (e) {
    return thunkAPI.rejectWithValue(e)
  }
})

export const getImageData = createAsyncThunk<
  { data: ImageData },
  { path: string; maxIndex?: number }
>(
  `${PLOT_DATA_SLICE_NAME}/getImageData`,
  async ({ path, maxIndex }, thunkAPI) => {
    try {
      const response = await axios.get(`${BASE_URL}/outputs/image/${path}`, {
        params: {
          max_index: maxIndex,
        },
      })
      return response.data
    } catch (e) {
      return thunkAPI.rejectWithValue(e)
    }
  },
)
