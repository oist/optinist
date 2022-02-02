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
}

export const DATA_TYPE_SET = {
  TIME_SERIES: 'timeSeries',
  HEAT_MAP: 'heatMap',
  IMAGE: 'image',
  CSV: 'csv',
  ROI: 'roi',
  SCATTER: 'scatter',
  BAR: 'bar',
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
  extends BaseDisplay<'timeSeries', TimeSeriesData> {}
export type TimeSeriesData = {
  [key: string]: {
    [key: number]: number
  }
}

export interface HeatMapDisplayData
  extends BaseDisplay<'heatMap', HeatMapData> {}
export type HeatMapData = number[][]

export interface ImageDisplayData extends BaseDisplay<'image', ImageData> {}
export type ImageData = number[][][]

export interface CsvDisplayData extends BaseDisplay<'csv', CsvData> {
  // columns: string[]
}
export type CsvData = number[][]

export interface RoiDisplayData extends BaseDisplay<'roi', RoiData> {}
export type RoiData = number[][][]

export interface ScatterDisplayData
  extends BaseDisplay<'scatter', ScatterData> {}
export type ScatterData = number[][][]

export interface BarDisplayData extends BaseDisplay<'bar', BarData> {}
export type BarData = number[][][]
