import { createSlice } from '@reduxjs/toolkit'
import { DisplayData, DISPLAY_DATA_SLICE_NAME } from './DisplayDataType'
import {
  getTableData,
  getHeatMapData,
  getImageData,
  getTimeSeriesData,
  getRoiData,
  getScatterData,
} from './DisplayDataActions'

const initialState: DisplayData = {
  timeSeries: {},
  heatMap: {},
  image: {},
  table: {},
  roi: {},
  scatter: {},
}

export const displayDataSlice = createSlice({
  name: DISPLAY_DATA_SLICE_NAME,
  initialState,
  reducers: {
    // deleteDisplayItem: (state, action: PayloadAction<number>) => {
    //   const itemId = action.payload
    //   delete state[itemId]
    //   if (itemId === state.selectedItemId) {
    //     state.selectedItemId = null
    //   }
    // },
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

export default displayDataSlice.reducer
