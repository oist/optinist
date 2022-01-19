import { RootState } from 'store/store'

const selectDisplayData = (state: RootState) => state.displayData

export const selectTimeSeriesData = (filePath: string) => (state: RootState) =>
  selectDisplayData(state).timeSeries[filePath].data

export const selectTimeSeriesPlotlyData =
  (filePath: string) => (state: RootState) => {
    let timeSeriesData = selectDisplayData(state).timeSeries[filePath].data
    if (timeSeriesData == null) {
      return []
    }
    return Object.keys(timeSeriesData['0']).map((_, i) => {
      return {
        name: `(${i})`,
        x: Object.keys(timeSeriesData),
        y: Object.values(timeSeriesData).map((value) => value[i]),
        visible: i === 0 ? true : 'legendonly',
      }
    })
  }

export const selectTimeSeriesDataIsInitialized =
  (filePath: string) => (state: RootState) =>
    Object.keys(selectDisplayData(state).timeSeries).includes(filePath)

export const selectTimeSeriesDataIsPending =
  (filePath: string) => (state: RootState) =>
    selectTimeSeriesDataIsInitialized(filePath)(state) &&
    selectDisplayData(state).timeSeries[filePath].pending

export const selectTimeSeriesDataIsFulfilled =
  (filePath: string) => (state: RootState) =>
    selectTimeSeriesDataIsInitialized(filePath)(state) &&
    selectDisplayData(state).timeSeries[filePath].fulfilled

export const selectTimeSeriesDataError =
  (filePath: string) => (state: RootState) =>
    selectTimeSeriesDataIsInitialized(filePath)(state)
      ? selectDisplayData(state).timeSeries[filePath].error
      : null

export const selectHeatMapData = (filePath: string) => (state: RootState) =>
  selectDisplayData(state).heatMap[filePath].data

export const selectHeatMapDataIsInitialized =
  (filePath: string) => (state: RootState) =>
    Object.keys(selectDisplayData(state).heatMap).includes(filePath)

export const selectHeatMapDataIsPending =
  (filePath: string) => (state: RootState) =>
    selectHeatMapDataIsInitialized(filePath)(state) &&
    selectDisplayData(state).heatMap[filePath].pending

export const selectHeatMapDataIsFulfilled =
  (filePath: string) => (state: RootState) =>
    selectHeatMapDataIsInitialized(filePath)(state) &&
    selectDisplayData(state).heatMap[filePath].fulfilled

export const selectHeatMapDataError =
  (filePath: string) => (state: RootState) =>
    selectHeatMapDataIsInitialized(filePath)(state)
      ? selectDisplayData(state).heatMap[filePath].error
      : null

export const selectImageData = (filePath: string) => (state: RootState) =>
  selectDisplayData(state).image[filePath]

export const selectImageDataIsInitialized =
  (filePath: string) => (state: RootState) =>
    Object.keys(selectDisplayData(state).image).includes(filePath)

export const selectImageDataError = (filePath: string) => (state: RootState) =>
  selectImageDataIsInitialized(filePath)(state)
    ? selectDisplayData(state).image[filePath].error
    : null

export const selectImageDataIsPending =
  (filePath: string) => (state: RootState) =>
    selectImageDataIsInitialized(filePath)(state) &&
    selectDisplayData(state).image[filePath].pending

export const selectImageDataIsFulfilled =
  (filePath: string) => (state: RootState) =>
    selectImageDataIsInitialized(filePath)(state) &&
    selectDisplayData(state).image[filePath].fulfilled

export const selectImageDataMaxIndex =
  (filePath: string) => (state: RootState) => {
    if (!selectImageDataIsPending(filePath)(state)) {
      return selectImageData(filePath)(state).data.length - 1
    } else {
      return undefined
    }
  }

export const selectActiveImageData =
  (filePath: string, activeIndex: number) => (state: RootState) => {
    return selectImageData(filePath)(state).data[activeIndex]
  }

export const selectTableData = (filePath: string) => (state: RootState) =>
  selectDisplayData(state).table[filePath].data

export const selectTableDataIsInitialized =
  (filePath: string) => (state: RootState) =>
    Object.keys(selectDisplayData(state).table).includes(filePath)

export const selectTableDataError = (filePath: string) => (state: RootState) =>
  selectTableDataIsInitialized(filePath)(state)
    ? selectDisplayData(state).table[filePath].error
    : null

export const selectTableDataIsPending =
  (filePath: string) => (state: RootState) =>
    selectTableDataIsInitialized(filePath)(state) &&
    selectDisplayData(state).table[filePath].pending

export const selectTableDataIsFulfilled =
  (filePath: string) => (state: RootState) =>
    selectTableDataIsInitialized(filePath)(state) &&
    selectDisplayData(state).table[filePath].fulfilled

export const selectTableDataColumns =
  (filePath: string) => (state: RootState) =>
    selectDisplayData(state).table[filePath].columns
