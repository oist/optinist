import { DATA_TYPE, DATA_TYPE_SET } from '../DisplayData/DisplayDataType'

export type VisualaizeItem = {
  selectedItemId: number | null
  items: {
    [itemId: number]: VisualaizeItemType
  }
}

export type VisualaizeItemType = DefaultSetItem | DisplayDataItem

export interface ItemBaseType<T extends VISUALIZE_ITEM_TYPE> {
  itemType: T
}

export type ColorType = {
  rgb: string
  offset: string
}

export const VISUALIZE_ITEM_TYPE_SET = {
  DEFAULT_SET: 'defaultSet',
  DISPLAY_DATA: 'displayData',
} as const

export type VISUALIZE_ITEM_TYPE =
  typeof VISUALIZE_ITEM_TYPE_SET[keyof typeof VISUALIZE_ITEM_TYPE_SET]

export type DisplayDataItem =
  | ImageItem
  | TimeSeriesItem
  | HeatMapItem
  | CsvItem
  | RoiItem
  | ScatterItem
  | BarItem
  | HDF5Item

export interface DisplayDataItemBaseType extends ItemBaseType<'displayData'> {
  filePath: string | null
  nodeId: string | null
  dataType: DATA_TYPE | null
}

export interface DefaultSetItem extends ItemBaseType<'defaultSet'> {
  imageItem: ImageItem
  timeSeriesItem: TimeSeriesItem
  heatMapItem: HeatMapItem
  otherItem?: DisplayDataItem
}

export interface ImageItem extends DisplayDataItemBaseType {
  dataType: typeof DATA_TYPE_SET.IMAGE
  activeIndex: number
  startIndex: number
  endIndex: number
  showticklabels: boolean
  zsmooth: string | boolean
  showline: boolean
  showgrid: boolean
  showscale: boolean
  colors: ColorType[]
  roiItem: RoiItem | null
}

export interface TimeSeriesItem extends DisplayDataItemBaseType {
  dataType: typeof DATA_TYPE_SET.TIME_SERIES
  offset: boolean
  span: number
  showgrid: boolean
  showline: boolean
  showticklabels: boolean
  zeroline: boolean
  xrange: {
    left: number | undefined
    right: number | undefined
  }
  maxIndex: number
}

export interface HeatMapItem extends DisplayDataItemBaseType {
  dataType: typeof DATA_TYPE_SET.HEAT_MAP
  showscale: boolean
  colors: ColorType[]
}
export interface CsvItem extends DisplayDataItemBaseType {
  dataType: typeof DATA_TYPE_SET.CSV
  setColumn: number | null
  setIndex: boolean
  transpose: boolean
}
export interface RoiItem extends DisplayDataItemBaseType {
  dataType: typeof DATA_TYPE_SET.ROI
  colors: ColorType[]
}

export interface ScatterItem extends DisplayDataItemBaseType {
  dataType: typeof DATA_TYPE_SET.SCATTER
  xIndex: number
  yIndex: number
}

export interface BarItem extends DisplayDataItemBaseType {
  dataType: typeof DATA_TYPE_SET.BAR
}

export interface HDF5Item extends DisplayDataItemBaseType {
  dataType: typeof DATA_TYPE_SET.HDF5
}
