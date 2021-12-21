import { createSlice, PayloadAction } from '@reduxjs/toolkit'
import { DisplayData, DISPLAY_DATA_SLICE_NAME } from './DisplayDataType'
import {
  getTableData,
  getHeatMapData,
  getImageData,
  getTimeSeriesData,
} from './DisplayDataActions'

const initialState: DisplayData = {
  timeSeries: {},
  heatMap: {},
  image: {},
  table: {},
}

export const displayDataSlice = createSlice({
  name: DISPLAY_DATA_SLICE_NAME,
  initialState,
  reducers: {
    incrementImageActiveIndex: (
      state,
      action: PayloadAction<{ path: string; amount?: number }>,
    ) => {
      const { path, amount } = action.payload
      const newIndex = state.image[path].activeIndex + (amount ?? 1)
      if (state.image[path].data[newIndex] != null) {
        state.image[path].activeIndex = newIndex
      }
    },
    decrementImageActiveIndex: (
      state,
      action: PayloadAction<{ path: string; amount?: number }>,
    ) => {
      const { path, amount } = action.payload
      const newIndex = state.image[path].activeIndex - (amount ?? 1)
      if (newIndex >= 0) {
        state.image[path].activeIndex = newIndex
      }
    },
  },
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
          activeIndex: 0,
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
          activeIndex: 0,
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
          activeIndex: 0,
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
  },
})

export const { decrementImageActiveIndex, incrementImageActiveIndex } =
  displayDataSlice.actions

export default displayDataSlice.reducer
