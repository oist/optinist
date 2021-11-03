import { createSlice } from '@reduxjs/toolkit'
import { toPlotDataKey } from './PlotDataUtils'
import { getHeatMapData, getTimeSeriesData } from './PlotDataAction'
import { PlotData, PLOT_DATA_SLICE_NAME } from './PlotDataType'

const initialState: PlotData = {
  timeSeriesDataMap: {},
  heatMapDataMap: {},
}

export const plotDataSlice = createSlice({
  name: PLOT_DATA_SLICE_NAME,
  initialState,
  reducers: {},
  extraReducers: (builder) => {
    builder
      .addCase(getTimeSeriesData.fulfilled, (state, action) => {
        const { nodeId, outputKey } = action.meta.arg
        const plotDataKey = toPlotDataKey(nodeId, outputKey)
        state.timeSeriesDataMap[plotDataKey] = action.payload
      })
      .addCase(getHeatMapData.fulfilled, (state, action) => {
        const { nodeId, outputKey } = action.meta.arg
        const plotDataKey = toPlotDataKey(nodeId, outputKey)
        state.heatMapDataMap[plotDataKey] = action.payload
      })
  },
})

// export const {} = plotDataSlice.actions

export default plotDataSlice.reducer
