import { createSlice, PayloadAction } from '@reduxjs/toolkit'
import { DisplayData, DISPLAY_DATA_SLICE_NAME } from './DisplayDataType'
import {
  getTableData,
  getHeatMapData,
  getImageData,
  getTimeSeriesData,
  getRoiData,
} from './DisplayDataActions'

const initialState: DisplayData = {
  timeSeries: {},
  heatMap: {},
  image: {},
  table: {},
  roi: {},
}

export const displayDataSlice = createSlice({
  name: DISPLAY_DATA_SLICE_NAME,
  initialState,
  reducers: {},
  extraReducers: (builder) => {
    builder
      .addCase(getTimeSeriesData.pending, (state, action) => {
        const { path } = action.meta.arg
        state.timeSeries[path] = {
          type: 'timeSeries',
          data: {},
          pending: true,
          fulfilled: false,
          error: null,
        }
      })
      .addCase(getTimeSeriesData.fulfilled, (state, action) => {
        const { path } = action.meta.arg
        state.timeSeries[path] = {
          type: 'timeSeries',
          data: action.payload.data,
          pending: false,
          fulfilled: true,
          error: null,
        }
      })
      .addCase(getTimeSeriesData.rejected, (state, action) => {
        const { path } = action.meta.arg
        state.timeSeries[path] = {
          type: 'timeSeries',
          data: {},
          pending: false,
          fulfilled: false,
          error: action.error.message ?? 'rejected',
        }
      })
      .addCase(getHeatMapData.pending, (state, action) => {
        const { path } = action.meta.arg
        state.heatMap[path] = {
          type: 'heatMap',
          data: [],
          pending: true,
          fulfilled: false,
          error: null,
        }
      })
      .addCase(getHeatMapData.fulfilled, (state, action) => {
        const { path } = action.meta.arg
        state.heatMap[path] = {
          type: 'heatMap',
          data: action.payload.data,
          pending: false,
          fulfilled: true,
          error: null,
        }
      })
      .addCase(getHeatMapData.rejected, (state, action) => {
        const { path } = action.meta.arg
        state.heatMap[path] = {
          type: 'heatMap',
          data: [],
          pending: false,
          fulfilled: false,
          error: action.error.message ?? 'rejected',
        }
      })
      .addCase(getImageData.pending, (state, action) => {
        const { path } = action.meta.arg
        state.image[path] = {
          type: 'image',
          // activeIndex: 0,
          data: [],
          pending: true,
          fulfilled: false,
          error: null,
        }
      })
      .addCase(getImageData.fulfilled, (state, action) => {
        const { path } = action.meta.arg
        state.image[path] = {
          type: 'image',
          // activeIndex: 0,
          data: action.payload.data,
          pending: false,
          fulfilled: true,
          error: null,
        }
      })
      .addCase(getImageData.rejected, (state, action) => {
        const { path } = action.meta.arg
        state.image[path] = {
          type: 'image',
          // activeIndex: 0,
          data: [],
          pending: false,
          fulfilled: false,
          error: action.error.message ?? 'rejected',
        }
      })
      .addCase(getTableData.pending, (state, action) => {
        const { path } = action.meta.arg
        state.table[path] = {
          type: 'table',
          columns: [],
          data: [],
          pending: true,
          fulfilled: false,
          error: null,
        }
      })
      .addCase(getTableData.fulfilled, (state, action) => {
        const { path } = action.meta.arg
        state.table[path] = {
          type: 'table',
          columns: action.payload.columns,
          data: action.payload.data,
          pending: false,
          fulfilled: true,
          error: null,
        }
      })
      .addCase(getTableData.rejected, (state, action) => {
        const { path } = action.meta.arg
        state.table[path] = {
          type: 'table',
          columns: [],
          data: [],
          pending: false,
          fulfilled: false,
          error: action.error.message ?? 'rejected',
        }
      })
      .addCase(getRoiData.pending, (state, action) => {
        const { path } = action.meta.arg
        state.roi[path] = {
          type: 'roi',
          data: [],
          pending: true,
          fulfilled: false,
          error: null,
        }
      })
      .addCase(getRoiData.fulfilled, (state, action) => {
        const { path } = action.meta.arg
        state.roi[path] = {
          type: 'roi',
          data: action.payload.data,
          pending: false,
          fulfilled: true,
          error: null,
        }
      })
      .addCase(getRoiData.rejected, (state, action) => {
        const { path } = action.meta.arg
        state.roi[path] = {
          type: 'roi',
          data: [],
          pending: false,
          fulfilled: false,
          error: action.error.message ?? 'rejected',
        }
      })
  },
})

export default displayDataSlice.reducer
