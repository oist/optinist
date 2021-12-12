import { createSlice, PayloadAction } from '@reduxjs/toolkit'
import {
  getHeatMapData,
  getTimeSeriesData,
  getImageData,
  getTableData,
} from './PlotDataAction'
import { PlotData, PLOT_DATA_SLICE_NAME } from './PlotDataType'

const initialState: PlotData = {
  timeSeriesDataMap: {},
  heatMapDataMap: {},
  imageDataMap: {},
  tableDataMap: {},
}

export const plotDataSlice = createSlice({
  name: PLOT_DATA_SLICE_NAME,
  initialState,
  reducers: {
    incrementImageActiveIndex: (
      state,
      action: PayloadAction<{ path: string; amount?: number }>,
    ) => {
      const { path, amount } = action.payload
      const newIndex = state.imageDataMap[path].activeIndex + (amount ?? 1)
      if (state.imageDataMap[path].data[newIndex] != null) {
        state.imageDataMap[path].activeIndex = newIndex
      }
    },
    decrementImageActiveIndex: (
      state,
      action: PayloadAction<{ path: string; amount?: number }>,
    ) => {
      const { path, amount } = action.payload
      const newIndex = state.imageDataMap[path].activeIndex - (amount ?? 1)
      if (newIndex >= 0) {
        state.imageDataMap[path].activeIndex = newIndex
      }
    },
  },
  extraReducers: (builder) => {
    builder
      .addCase(getTimeSeriesData.fulfilled, (state, action) => {
        const { path } = action.meta.arg
        state.timeSeriesDataMap[path] = {
          data: action.payload.data,
        }
      })
      .addCase(getHeatMapData.fulfilled, (state, action) => {
        const { path } = action.meta.arg
        state.heatMapDataMap[path] = {
          data: action.payload.data,
        }
      })
      .addCase(getImageData.fulfilled, (state, action) => {
        const { path } = action.meta.arg
        state.imageDataMap[path] = {
          activeIndex: 0,
          data: action.payload.data,
        }
      })
      .addCase(getTableData.fulfilled, (state, action) => {
        const { path } = action.meta.arg
        state.tableDataMap[path] = {
          columns: action.payload.columns,
          data: action.payload.data,
        }
      })
  },
})

export const { incrementImageActiveIndex, decrementImageActiveIndex } =
  plotDataSlice.actions

export default plotDataSlice.reducer
