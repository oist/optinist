import { DATA_TYPE, DATA_TYPE_SET } from '../DisplayData/DisplayDataType'

export const VISUALIZE_ITEM_SLICE_NAME = 'visualaizeItem'

export type VisualaizeItem = {
  selectedItemId: number | null
  items: {
    [itemId: number]: VisualaizeItemType
  }
  layout: ItemLayout
}

export type ItemLayout = number[][] // itemIdをrow,columnで並べる

export type VisualaizeItemType = DisplayDataItem

export interface ItemBaseType<T extends VISUALIZE_ITEM_TYPE> {
  itemType: T
  width: number
  height: number
}

export type ColorType = {
  rgb: string
  offset: string
}

export const VISUALIZE_ITEM_TYPE_SET = {
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
  | HTMLItem
  | FluoItem
  | BehaviorItem

export interface DisplayDataItemBaseType extends ItemBaseType<'displayData'> {
  filePath: string | null
  nodeId: string | null
  dataType: DATA_TYPE | null
  isWorkflowDialog: boolean
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
  roiAlpha: number
  duration: number
  clickedDataId?: number
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
  displayNumbers: number[]
  refImageItemId?: number
  checkedList: boolean[]
}

export interface HeatMapItem extends DisplayDataItemBaseType {
  dataType: typeof DATA_TYPE_SET.HEAT_MAP
  showscale: boolean
  colors: ColorType[]
}
export interface CsvItem extends DisplayDataItemBaseType {
  dataType: typeof DATA_TYPE_SET.CSV
  setHeader: number | null
  setIndex: boolean
  transpose: boolean
}
export interface RoiItem extends DisplayDataItemBaseType {
  dataType: typeof DATA_TYPE_SET.ROI
  // colors: ColorType[]
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

export interface HTMLItem extends DisplayDataItemBaseType {
  dataType: typeof DATA_TYPE_SET.HTML
}

export interface FluoItem extends DisplayDataItemBaseType {
  dataType: typeof DATA_TYPE_SET.FLUO
}

export interface BehaviorItem extends DisplayDataItemBaseType {
  dataType: typeof DATA_TYPE_SET.BEHAVIOR
}
