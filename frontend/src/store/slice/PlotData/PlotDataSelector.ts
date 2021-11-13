import { RootState } from '../../store'

export const plotDataSelector = (state: RootState) => state.plotData

export const timeSeriesDataSelector = (path: string) => (state: RootState) => {
  if (Object.keys(state.plotData.timeSeriesDataMap).includes(path)) {
    return state.plotData.timeSeriesDataMap[path].data
  } else {
    return undefined
  }
}

export const timeSeriesDataIsLoadedSelector =
  (path: string) => (state: RootState) => {
    return Object.keys(state.plotData.timeSeriesDataMap).includes(path)
  }

export const heatMapDataSelector = (path: string) => (state: RootState) => {
  if (Object.keys(state.plotData.heatMapDataMap).includes(path)) {
    return state.plotData.heatMapDataMap[path].data
  } else {
    return undefined
  }
}

export const heatMapDataIsLoadedSelector =
  (path: string) => (state: RootState) => {
    return Object.keys(state.plotData.heatMapDataMap).includes(path)
  }

export const activeImageDataSelector = (path: string) => (state: RootState) => {
  if (Object.keys(state.plotData.imageDataMap).includes(path)) {
    const index = state.plotData.imageDataMap[path].activeIndex
    return state.plotData.imageDataMap[path].data[index]
  } else {
    return undefined
  }
}

export const imageDataIsLoadedSelector =
  (path: string) => (state: RootState) => {
    return Object.keys(state.plotData.imageDataMap).includes(path)
  }

export const imageDataActiveIndexselector =
  (path: string) => (state: RootState) => {
    if (Object.keys(state.plotData.imageDataMap).includes(path)) {
      return state.plotData.imageDataMap[path].activeIndex
    } else {
      return undefined
    }
  }

export const imageDataMaxIndexSelector =
  (path: string) => (state: RootState) => {
    if (Object.keys(state.plotData.imageDataMap).includes(path)) {
      return state.plotData.imageDataMap[path].data.length - 1
    } else {
      return undefined
    }
  }
