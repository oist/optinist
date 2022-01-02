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
  table: {
    [filePath: string]: TableDisplayData
  }
}

export const DATA_TYPE_SET = {
  TIME_SERIES: 'timeSeries',
  HEAT_MAP: 'heatMap',
  IMAGE: 'image',
  TABLE: 'table',
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

export interface ImageDisplayData extends BaseDisplay<'image', ImageData> {
  activeIndex: number
}
export type ImageData = number[][][]

export interface TableDisplayData extends BaseDisplay<'table', TableData> {
  columns: string[]
}
export type TableData = number[][]
