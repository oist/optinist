import { createAsyncThunk } from '@reduxjs/toolkit'
import { PLOT_DATA_SLICE_NAME } from './PlotDataType'
import axios from 'axios'

import { TimeSeriesData } from './PlotDataType'
import { HeatMapData } from './PlotDataType'
import { BASE_URL } from 'const/API'

export const getTimeSeriesData = createAsyncThunk<
  TimeSeriesData,
  { nodeId: string; outputKey: string; path: string }
>(`${PLOT_DATA_SLICE_NAME}/getTimeSeriesData`, async ({ path }, thunkAPI) => {
  try {
    const response = await axios.get(`${BASE_URL}/api/outputs/${path}`)
    return response.data
  } catch (e) {
    return thunkAPI.rejectWithValue(e)
  }
})

export const getHeatMapData = createAsyncThunk<
  HeatMapData,
  { nodeId: string; outputKey: string; path: string }
>(`${PLOT_DATA_SLICE_NAME}/getHeatMapData`, async ({ path }, thunkAPI) => {
  try {
    const response = await axios.get(`${BASE_URL}/api/outputs/${path}`)
    return response.data
  } catch (e) {
    return thunkAPI.rejectWithValue(e)
  }
})
