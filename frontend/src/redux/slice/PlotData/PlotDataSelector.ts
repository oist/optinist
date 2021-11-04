import { RootState } from '../../store'
import { PlotDataKey } from './PlotDataType'

export const plotDataSelector = (state: RootState) => state.plotData

export const timeSeriesDataByKeySelector =
  (id: PlotDataKey) => (state: RootState) => {
    if (Object.keys(state.plotData.timeSeriesDataMap).includes(id)) {
      return state.plotData.timeSeriesDataMap[id]
    } else {
      return undefined
    }
  }

export const timeSeriesDataIsLoadedByKeySelector =
  (key: PlotDataKey) => (state: RootState) => {
    return Object.keys(state.plotData.timeSeriesDataMap).includes(key)
  }

export const heatMapDataByKeySelector =
  (key: PlotDataKey) => (state: RootState) => {
    if (Object.keys(state.plotData.heatMapDataMap).includes(key)) {
      return state.plotData.heatMapDataMap[key]
    } else {
      return undefined
    }
  }

export const heatMapDataIsLoadedByKeySelector =
  (key: PlotDataKey) => (state: RootState) => {
    return Object.keys(state.plotData.heatMapDataMap).includes(key)
  }
