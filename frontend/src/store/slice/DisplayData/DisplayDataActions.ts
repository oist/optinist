import { createAsyncThunk } from '@reduxjs/toolkit'
import {
  TimeSeriesData,
  CsvData,
  RoiData,
  ScatterData,
  BarData,
  HTMLData,
  ImageData,
  HeatMapData,
  getTimeSeriesInitDataApi,
  getTimeSeriesDataByIdApi,
  getTimeSeriesAllDataApi,
  getHeatMapDataApi,
  getImageDataApi,
  getCsvDataApi,
  getRoiDataApi,
  getScatterDataApi,
  getBarDataApi,
  getHTMLDataApi,
} from 'api/outputs/Outputs'
import { DISPLAY_DATA_SLICE_NAME } from './DisplayDataType'

export const getTimeSeriesInitData = createAsyncThunk<
  { data: TimeSeriesData; xrange: number[]; std: TimeSeriesData },
  { path: string; itemId: number }
>(
  `${DISPLAY_DATA_SLICE_NAME}/getTimeSeriesInitData`,
  async ({ path }, thunkAPI) => {
    try {
      const response = await getTimeSeriesInitDataApi(path)
      return response
    } catch (e) {
      return thunkAPI.rejectWithValue(e)
    }
  },
)

export const getTimeSeriesDataById = createAsyncThunk<
  { data: TimeSeriesData; xrange: number[]; std: TimeSeriesData },
  { path: string; index: string }
>(
  `${DISPLAY_DATA_SLICE_NAME}/getTimeSeriesDataById`,
  async ({ path, index }, thunkAPI) => {
    try {
      const response = await getTimeSeriesDataByIdApi(path, index)
      return response
    } catch (e) {
      return thunkAPI.rejectWithValue(e)
    }
  },
)

export const getTimeSeriesAllData = createAsyncThunk<
  { data: TimeSeriesData; xrange: number[]; std: TimeSeriesData },
  { path: string }
>(
  `${DISPLAY_DATA_SLICE_NAME}/getTimeSeriesAllData`,
  async ({ path }, thunkAPI) => {
    try {
      const response = await getTimeSeriesAllDataApi(path)
      return response
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
    const response = await getHeatMapDataApi(path)
    return response
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
      const response = await getImageDataApi(path, { startIndex, endIndex })
      return response
    } catch (e) {
      return thunkAPI.rejectWithValue(e)
    }
  },
)

export const getCsvData = createAsyncThunk<
  {
    data: CsvData
  },
  { path: string }
>(`${DISPLAY_DATA_SLICE_NAME}/getCsvData`, async ({ path }, thunkAPI) => {
  try {
    const response = await getCsvDataApi(path)
    return response
  } catch (e) {
    return thunkAPI.rejectWithValue(e)
  }
})

export const getRoiData = createAsyncThunk<{ data: RoiData }, { path: string }>(
  `${DISPLAY_DATA_SLICE_NAME}/getRoiData`,
  async ({ path }, thunkAPI) => {
    try {
      const response = await getRoiDataApi(path)
      return response
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
    const response = await getScatterDataApi(path)
    return response
  } catch (e) {
    return thunkAPI.rejectWithValue(e)
  }
})

export const getBarData = createAsyncThunk<{ data: BarData }, { path: string }>(
  `${DISPLAY_DATA_SLICE_NAME}/getBarData`,
  async ({ path }, thunkAPI) => {
    try {
      const response = await getBarDataApi(path)
      return response
    } catch (e) {
      return thunkAPI.rejectWithValue(e)
    }
  },
)

export const getHTMLData = createAsyncThunk<
  { data: HTMLData },
  { path: string }
>(`${DISPLAY_DATA_SLICE_NAME}/getHTMLData`, async ({ path }, thunkAPI) => {
  try {
    const response = await getHTMLDataApi(path)
    return response
  } catch (e) {
    return thunkAPI.rejectWithValue(e)
  }
})
