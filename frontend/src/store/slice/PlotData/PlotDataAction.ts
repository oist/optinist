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
    const response = await axios.get(`${BASE_URL}/api/outputs/${path}`)
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
    const response = await axios.get(`${BASE_URL}/api/outputs/${path}`)
    return response.data
  } catch (e) {
    return thunkAPI.rejectWithValue(e)
  }
})

export const getImageData = createAsyncThunk<
  { data: ImageData },
  { path: string }
>(`${PLOT_DATA_SLICE_NAME}/getImageData`, async ({ path }, thunkAPI) => {
  try {
    const response = await axios.get(`${BASE_URL}/api/outputs/${path}`)
    return response.data
  } catch (e) {
    return thunkAPI.rejectWithValue(e)
  }
})
