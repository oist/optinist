import {
  HeatMapData,
  TimeSeriesData,
  BarData,
  CsvData,
  RoiData,
  ImageData,
  ScatterData,
  HTMLData,
} from 'api/outputs/Outputs'

export const DISPLAY_DATA_SLICE_NAME = 'displayData'

export type DisplayData = {
  timeSeries: {
    [filePath: string]: TimeSeriesDisplayData
  }
  heatMap: {
    [filePath: string]: HeatMapDisplayData
  }
  image: {
    [filePath: string]: ImageDisplayData
  }
  csv: {
    [filePath: string]: CsvDisplayData
  }
  roi: {
    [filePath: string]: RoiDisplayData
  }
  scatter: {
    [filePath: string]: ScatterDisplayData
  }
  bar: {
    [filePath: string]: BarDisplayData
  }
  // hdf5: {
  //   [filePath: string]: HDF5DisplayData
  // }
  html: {
    [filePath: string]: HTMLDisplayData
  }
}

export const DATA_TYPE_SET = {
  TIME_SERIES: 'timeSeries',
  HEAT_MAP: 'heatMap',
  IMAGE: 'image',
  CSV: 'csv',
  ROI: 'roi',
  SCATTER: 'scatter',
  BAR: 'bar',
  HDF5: 'hdf5',
  HTML: 'html',
  FLUO: 'fluo',
  BEHAVIOR: 'behavior',
} as const

export type DATA_TYPE = typeof DATA_TYPE_SET[keyof typeof DATA_TYPE_SET]

interface BaseDisplay<T extends DATA_TYPE, Data> {
  type: T
  data: Data
  pending: boolean
  error: string | null
  fulfilled: boolean
}

export interface TimeSeriesDisplayData
  extends BaseDisplay<'timeSeries', TimeSeriesData> {
  xrange: number[]
  std: TimeSeriesData
}

export interface HeatMapDisplayData
  extends BaseDisplay<'heatMap', HeatMapData> {
  columns: string[]
  index: string[]
}

export interface ImageDisplayData extends BaseDisplay<'image', ImageData> {}

export interface CsvDisplayData extends BaseDisplay<'csv', CsvData> {
  // columns: string[]
}

export interface RoiDisplayData extends BaseDisplay<'roi', RoiData> {
  roiUniqueList: string[]
}

export interface ScatterDisplayData
  extends BaseDisplay<'scatter', ScatterData> {}

export interface BarDisplayData extends BaseDisplay<'bar', BarData> {
  columns: string[]
  index: string[]
}

export interface HTMLDisplayData extends BaseDisplay<'html', HTMLData> {}
