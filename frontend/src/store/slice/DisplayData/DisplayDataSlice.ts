import { useSelector } from 'react-redux'
import { createSlice, PayloadAction } from '@reduxjs/toolkit'
import { DisplayData, DISPLAY_DATA_SLICE_NAME } from './DisplayDataType'
import {
  getCsvData,
  getHeatMapData,
  getImageData,
  getTimeSeriesData,
  getRoiData,
  getScatterData,
} from './DisplayDataActions'
import {
  selectVisualizeDataFilePath,
  selectVisualizeDataType,
} from '../VisualizeItem/VisualizeItemSelectors'

const initialState: DisplayData = {
  timeSeries: {},
  heatMap: {},
  image: {},
  csv: {},
  roi: {},
  scatter: {},
}

export const displayDataSlice = createSlice({
  name: DISPLAY_DATA_SLICE_NAME,
  initialState,
  reducers: {
    deleteDisplayImageItem: (
      state,
      action: PayloadAction<{
        filePath: string
      }>,
    ) => {
      delete state.image[action.payload.filePath]
    },
    deleteDisplayTimeSeriesItem: (
      state,
      action: PayloadAction<{
        filePath: string
      }>,
    ) => {
      delete state.timeSeries[action.payload.filePath]
    },
    deleteDisplayCsvItem: (
      state,
      action: PayloadAction<{
        filePath: string
      }>,
    ) => {
      delete state.csv[action.payload.filePath]
    },
    deleteDisplayHeatMapItem: (
      state,
      action: PayloadAction<{
        filePath: string
      }>,
    ) => {
      delete state.heatMap[action.payload.filePath]
    },
    deleteDisplayRoiItem: (
      state,
      action: PayloadAction<{
        filePath: string
      }>,
    ) => {
      delete state.roi[action.payload.filePath]
    },
    deleteDisplayScatterItem: (
      state,
      action: PayloadAction<{
        filePath: string
      }>,
    ) => {
      delete state.scatter[action.payload.filePath]
    },
  },
  extraReducers: (builder) => {
    builder
      .addCase(getTimeSeriesData.pending, (state, action) => {
        const { path } = action.meta.arg
        if (!state.timeSeries.hasOwnProperty(path)) {
          state.timeSeries[path] = {
            type: 'timeSeries',
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
      .addCase(getTimeSeriesData.fulfilled, (state, action) => {
        const { path, index } = action.meta.arg
        state.timeSeries[path].pending = false
        state.timeSeries[path].fulfilled = true
        state.timeSeries[path].error = null
        if (Object.keys(state.timeSeries[path].data).length === 0) {
          state.timeSeries[path].data = action.payload.data
        } else {
          state.timeSeries[path].data[index] = action.payload.data[index]
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
          // columns: [],
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
          // columns: action.payload.columns,
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
          // columns: [],
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
  },
})

export const {
  deleteDisplayImageItem,
  deleteDisplayTimeSeriesItem,
  deleteDisplayCsvItem,
  deleteDisplayHeatMapItem,
  deleteDisplayRoiItem,
  deleteDisplayScatterItem,
} = displayDataSlice.actions

export default displayDataSlice.reducer
