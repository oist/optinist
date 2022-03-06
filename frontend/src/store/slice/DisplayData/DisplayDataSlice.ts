import { createSlice, PayloadAction } from '@reduxjs/toolkit'
import {
  DATA_TYPE,
  DATA_TYPE_SET,
  DisplayData,
  DISPLAY_DATA_SLICE_NAME,
} from './DisplayDataType'
import {
  getCsvData,
  getHeatMapData,
  getImageData,
  getTimeSeriesDataById,
  getTimeSeriesAllData,
  getRoiData,
  getScatterData,
  getBarData,
  getHTMLData,
} from './DisplayDataActions'

const initialState: DisplayData = {
  timeSeries: {},
  heatMap: {},
  image: {},
  csv: {},
  roi: {},
  scatter: {},
  bar: {},
  hdf5: {},
  html: {},
}

export const displayDataSlice = createSlice({
  name: DISPLAY_DATA_SLICE_NAME,
  initialState,
  reducers: {
    deleteDisplayItem: (
      state,
      action: PayloadAction<{
        dataType: DATA_TYPE
        filePath: string | null
      }>,
    ) => {
      const { dataType, filePath } = action.payload
      if (filePath !== null) {
        if (dataType === DATA_TYPE_SET.IMAGE) {
          delete state.image[filePath]
        } else if (dataType === DATA_TYPE_SET.TIME_SERIES) {
          delete state.timeSeries[filePath]
        } else if (dataType === DATA_TYPE_SET.CSV) {
          delete state.csv[filePath]
        } else if (dataType === DATA_TYPE_SET.HEAT_MAP) {
          delete state.heatMap[filePath]
        } else if (dataType === DATA_TYPE_SET.ROI) {
          delete state.roi[filePath]
        } else if (dataType === DATA_TYPE_SET.SCATTER) {
          delete state.scatter[filePath]
        } else if (dataType === DATA_TYPE_SET.BAR) {
          delete state.bar[filePath]
        } else if (dataType === DATA_TYPE_SET.HTML) {
          delete state.html[filePath]
        }
      }
    },
  },
  extraReducers: (builder) => {
    builder
      .addCase(getTimeSeriesDataById.pending, (state, action) => {
        const { path } = action.meta.arg
        if (!state.timeSeries.hasOwnProperty(path)) {
          state.timeSeries[path] = {
            type: 'timeSeries',
            xrange: [],
            data: {},
            pending: true,
            fulfilled: false,
            error: null,
          }
        } else {
          state.timeSeries[path].pending = true
          state.timeSeries[path].fulfilled = false
          state.timeSeries[path].error = null
        }
      })
      .addCase(getTimeSeriesDataById.fulfilled, (state, action) => {
        const { path, index } = action.meta.arg
        state.timeSeries[path].pending = false
        state.timeSeries[path].fulfilled = true
        state.timeSeries[path].error = null
        if (Object.keys(state.timeSeries[path].data).length === 0) {
          state.timeSeries[path].xrange = action.payload.data.xrange
          state.timeSeries[path].data = action.payload.data.data
        } else {
          state.timeSeries[path].data[index] = action.payload.data.data[index]
        }
      })
      .addCase(getTimeSeriesDataById.rejected, (state, action) => {
        const { path } = action.meta.arg
        state.timeSeries[path] = {
          type: 'timeSeries',
          xrange: [],
          data: {},
          pending: false,
          fulfilled: false,
          error: action.error.message ?? 'rejected',
        }
      })
      .addCase(getTimeSeriesAllData.pending, (state, action) => {
        const { path } = action.meta.arg
        if (!state.timeSeries.hasOwnProperty(path)) {
          state.timeSeries[path] = {
            type: 'timeSeries',
            xrange: [],
            data: {},
            pending: true,
            fulfilled: false,
            error: null,
          }
        } else {
          state.timeSeries[path].pending = true
          state.timeSeries[path].fulfilled = false
          state.timeSeries[path].error = null
        }
      })
      .addCase(getTimeSeriesAllData.rejected, (state, action) => {
        const { path } = action.meta.arg
        state.timeSeries[path] = {
          type: 'timeSeries',
          xrange: [],
          data: {},
          pending: false,
          fulfilled: false,
          error: action.error.message ?? 'rejected',
        }
      })
      .addCase(getTimeSeriesAllData.fulfilled, (state, action) => {
        const { path } = action.meta.arg
        state.timeSeries[path].pending = false
        state.timeSeries[path].fulfilled = true
        state.timeSeries[path].error = null
        state.timeSeries[path].xrange = action.payload.data.xrange
        state.timeSeries[path].data = action.payload.data.data
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
          data: [],
          pending: false,
          fulfilled: false,
          error: action.error.message ?? 'rejected',
        }
      })
      .addCase(getCsvData.pending, (state, action) => {
        const { path } = action.meta.arg
        state.csv[path] = {
          type: 'csv',
          data: [],
          pending: true,
          fulfilled: false,
          error: null,
        }
      })
      .addCase(getCsvData.fulfilled, (state, action) => {
        const { path } = action.meta.arg
        state.csv[path] = {
          type: 'csv',
          data: action.payload.data,
          pending: false,
          fulfilled: true,
          error: null,
        }
      })
      .addCase(getCsvData.rejected, (state, action) => {
        const { path } = action.meta.arg
        state.csv[path] = {
          type: 'csv',
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
      .addCase(getScatterData.pending, (state, action) => {
        const { path } = action.meta.arg
        state.scatter[path] = {
          type: 'scatter',
          data: [],
          pending: true,
          fulfilled: false,
          error: null,
        }
      })
      .addCase(getScatterData.fulfilled, (state, action) => {
        const { path } = action.meta.arg
        state.scatter[path] = {
          type: 'scatter',
          data: action.payload.data,
          pending: false,
          fulfilled: true,
          error: null,
        }
      })
      .addCase(getScatterData.rejected, (state, action) => {
        const { path } = action.meta.arg
        state.scatter[path] = {
          type: 'scatter',
          data: [],
          pending: false,
          fulfilled: false,
          error: action.error.message ?? 'rejected',
        }
      })
      .addCase(getBarData.pending, (state, action) => {
        const { path } = action.meta.arg
        state.bar[path] = {
          type: 'bar',
          data: [],
          pending: true,
          fulfilled: false,
          error: null,
        }
      })
      .addCase(getBarData.fulfilled, (state, action) => {
        const { path } = action.meta.arg
        state.bar[path] = {
          type: 'bar',
          data: action.payload.data,
          pending: false,
          fulfilled: true,
          error: null,
        }
      })
      .addCase(getBarData.rejected, (state, action) => {
        const { path } = action.meta.arg
        state.bar[path] = {
          type: 'bar',
          data: [],
          pending: false,
          fulfilled: false,
          error: action.error.message ?? 'rejected',
        }
      })
      .addCase(getHTMLData.pending, (state, action) => {
        const { path } = action.meta.arg
        state.html[path] = {
          type: 'html',
          data: '',
          pending: true,
          fulfilled: false,
          error: null,
        }
      })
      .addCase(getHTMLData.fulfilled, (state, action) => {
        const { path } = action.meta.arg
        state.html[path] = {
          type: 'html',
          data: action.payload.data,
          pending: false,
          fulfilled: true,
          error: null,
        }
      })
      .addCase(getHTMLData.rejected, (state, action) => {
        const { path } = action.meta.arg
        state.html[path] = {
          type: 'html',
          data: '',
          pending: false,
          fulfilled: false,
          error: action.error.message ?? 'rejected',
        }
      })
  },
})

export const { deleteDisplayItem } = displayDataSlice.actions

export default displayDataSlice.reducer
