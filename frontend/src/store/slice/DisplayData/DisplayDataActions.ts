import { createAsyncThunk } from "@reduxjs/toolkit"

import {
  TimeSeriesData,
  CsvData,
  RoiData,
  ScatterData,
  BarData,
  HTMLData,
  ImageData,
  HeatMapData,
  HistogramData,
  LineData,
  PieData,
  PolarData,
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
  getHistogramDataApi,
  getLineDataApi,
  getPieDataApi,
  getPolarDataApi,
  cancelRoiApi,
  addRoiApi,
  mergeRoiApi,
  deleteRoiApi,
  commitRoiApi,
  getStatusRoi,
} from "api/outputs/Outputs"
import { StatusROI } from "components/Workspace/Visualize/Plot/ImagePlot"
import {
  PlotMetaData,
  DISPLAY_DATA_SLICE_NAME,
} from "store/slice/DisplayData/DisplayDataType"

export const getTimeSeriesInitData = createAsyncThunk<
  {
    data: TimeSeriesData
    xrange: string[]
    std: TimeSeriesData
    meta?: PlotMetaData
  },
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
  {
    data: TimeSeriesData
    xrange: string[]
    std: TimeSeriesData
    meta?: PlotMetaData
  },
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
  {
    data: TimeSeriesData
    xrange: string[]
    std: TimeSeriesData
    meta?: PlotMetaData
  },
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
  {
    data: HeatMapData
    columns: string[]
    index: string[]
    meta?: PlotMetaData
  },
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
  { data: ImageData; meta?: PlotMetaData },
  { path: string; workspaceId: number; startIndex?: number; endIndex?: number }
>(
  `${DISPLAY_DATA_SLICE_NAME}/getImageData`,
  async ({ path, startIndex, endIndex, workspaceId }, thunkAPI) => {
    try {
      const response = await getImageDataApi(path, {
        workspaceId,
        startIndex,
        endIndex,
      })
      return response
    } catch (e) {
      return thunkAPI.rejectWithValue(e)
    }
  },
)

export const getCsvData = createAsyncThunk<
  { data: CsvData; meta?: PlotMetaData },
  { path: string; workspaceId: number }
>(
  `${DISPLAY_DATA_SLICE_NAME}/getCsvData`,
  async ({ path, workspaceId }, thunkAPI) => {
    try {
      const response = await getCsvDataApi(path, { workspaceId })
      return response
    } catch (e) {
      return thunkAPI.rejectWithValue(e)
    }
  },
)

export const getRoiData = createAsyncThunk<
  { data: RoiData; meta?: PlotMetaData },
  { path: string; workspaceId: number }
>(
  `${DISPLAY_DATA_SLICE_NAME}/getRoiData`,
  async ({ path, workspaceId }, thunkAPI) => {
    try {
      const response = await getRoiDataApi(path, { workspaceId })
      return response
    } catch (e) {
      return thunkAPI.rejectWithValue(e)
    }
  },
)

export const cancelRoi = createAsyncThunk<
  { data: HTMLData; meta?: PlotMetaData },
  { path: string | string[]; workspaceId: number }
>(
  `${DISPLAY_DATA_SLICE_NAME}/cancelRoi`,
  async ({ path, workspaceId }, thunkAPI) => {
    try {
      const response = await cancelRoiApi(path, workspaceId)
      return response
    } catch (e) {
      return thunkAPI.rejectWithValue(e)
    }
  },
)

export const addRoi = createAsyncThunk<
  { data: HTMLData; meta?: PlotMetaData },
  {
    path: string
    workspaceId: number
    data: { posx: number; posy: number; sizex: number; sizey: number }
  }
>(
  `${DISPLAY_DATA_SLICE_NAME}/addRoi`,
  async ({ path, workspaceId, data }, thunkAPI) => {
    const { dispatch } = thunkAPI
    try {
      const response = await addRoiApi(path, workspaceId, data)
      await dispatch(getStatus({ path, workspaceId }))
      await dispatch(getRoiData({ path, workspaceId }))
      return response
    } catch (e) {
      return thunkAPI.rejectWithValue(e)
    }
  },
)

export const mergeRoi = createAsyncThunk<
  { data: HTMLData; meta?: PlotMetaData },
  { path: string; workspaceId: number; data: { ids: number[] } }
>(
  `${DISPLAY_DATA_SLICE_NAME}/mergeROi`,
  async ({ path, workspaceId, data }, thunkAPI) => {
    const { dispatch } = thunkAPI
    try {
      const response = await mergeRoiApi(path, workspaceId, data)
      dispatch(getStatus({ path, workspaceId }))
      return response
    } catch (e) {
      return thunkAPI.rejectWithValue(e)
    }
  },
)

export const deleteRoi = createAsyncThunk<
  { data: HTMLData; meta?: PlotMetaData },
  { path: string; workspaceId: number; data: { ids: number[] } }
>(
  `${DISPLAY_DATA_SLICE_NAME}/deleteRoi`,
  async ({ path, workspaceId, data }, thunkAPI) => {
    const { dispatch } = thunkAPI
    try {
      const response = await deleteRoiApi(path, workspaceId, data)
      dispatch(getStatus({ path, workspaceId }))
      return response
    } catch (e) {
      return thunkAPI.rejectWithValue(e)
    }
  },
)

export const commitRoi = createAsyncThunk<
  boolean,
  { path: string; workspaceId: number }
>(
  `${DISPLAY_DATA_SLICE_NAME}/commitRoi`,
  async ({ path, workspaceId }, thunkAPI) => {
    const { dispatch } = thunkAPI
    try {
      const response = await commitRoiApi(path, workspaceId)
      dispatch(getStatus({ path, workspaceId }))
      return response
    } catch (e) {
      return thunkAPI.rejectWithValue(e)
    }
  },
)

export const getStatus = createAsyncThunk<
  StatusROI,
  { path: string; workspaceId: number }
>(
  `${DISPLAY_DATA_SLICE_NAME}/getStatusRoi`,
  async ({ path, workspaceId }, thunkAPI) => {
    try {
      const response = await getStatusRoi(path, workspaceId)
      return response
    } catch (e) {
      return thunkAPI.rejectWithValue(e)
    }
  },
)

export const getScatterData = createAsyncThunk<
  { data: ScatterData; meta?: PlotMetaData },
  { path: string }
>(`${DISPLAY_DATA_SLICE_NAME}/getScatterData`, async ({ path }, thunkAPI) => {
  try {
    const response = await getScatterDataApi(path)
    return response
  } catch (e) {
    return thunkAPI.rejectWithValue(e)
  }
})

export const getBarData = createAsyncThunk<
  { data: BarData; columns: string[]; index: string[]; meta?: PlotMetaData },
  { path: string }
>(`${DISPLAY_DATA_SLICE_NAME}/getBarData`, async ({ path }, thunkAPI) => {
  try {
    const response = await getBarDataApi(path)
    return response
  } catch (e) {
    return thunkAPI.rejectWithValue(e)
  }
})

export const getHTMLData = createAsyncThunk<
  { data: HTMLData; meta?: PlotMetaData },
  { path: string }
>(`${DISPLAY_DATA_SLICE_NAME}/getHTMLData`, async ({ path }, thunkAPI) => {
  try {
    const response = await getHTMLDataApi(path)
    return response
  } catch (e) {
    return thunkAPI.rejectWithValue(e)
  }
})

export const getHistogramData = createAsyncThunk<
  { data: HistogramData; meta?: PlotMetaData },
  { path: string }
>(`${DISPLAY_DATA_SLICE_NAME}/getHistogramData`, async ({ path }, thunkAPI) => {
  try {
    const response = await getHistogramDataApi(path)
    return response
  } catch (e) {
    return thunkAPI.rejectWithValue(e)
  }
})

export const getLineData = createAsyncThunk<
  { data: LineData; columns: number[]; index: number[]; meta?: PlotMetaData },
  { path: string }
>(`${DISPLAY_DATA_SLICE_NAME}/getLineData`, async ({ path }, thunkAPI) => {
  try {
    const response = await getLineDataApi(path)
    return response
  } catch (e) {
    return thunkAPI.rejectWithValue(e)
  }
})

export const getPieData = createAsyncThunk<
  { data: PieData; columns: string[]; meta?: PlotMetaData },
  { path: string }
>(`${DISPLAY_DATA_SLICE_NAME}/getPieData`, async ({ path }, thunkAPI) => {
  try {
    const response = await getPieDataApi(path)
    return response
  } catch (e) {
    return thunkAPI.rejectWithValue(e)
  }
})

export const getPolarData = createAsyncThunk<
  { data: PolarData; columns: number[]; index: number[]; meta?: PlotMetaData },
  { path: string }
>(`${DISPLAY_DATA_SLICE_NAME}/getPolarData`, async ({ path }, thunkAPI) => {
  try {
    const response = await getPolarDataApi(path)
    return response
  } catch (e) {
    return thunkAPI.rejectWithValue(e)
  }
})
