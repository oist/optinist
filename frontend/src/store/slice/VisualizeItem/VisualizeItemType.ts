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

export const VISUALIZE_ITEM_TYPE_SET = {
  DEFAULT_SET: 'defaultSet',
  DISPLAY_DATA: 'displayData',
} as const

export type VISUALIZE_ITEM_TYPE =
  typeof VISUALIZE_ITEM_TYPE_SET[keyof typeof VISUALIZE_ITEM_TYPE_SET]

export interface DisplayDataItem extends ItemBaseType<'displayData'> {
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

interface ImageItem extends DisplayDataItem {
  dataType: typeof DATA_TYPE_SET.IMAGE
}
interface TimeSeriesItem extends DisplayDataItem {
  dataType: typeof DATA_TYPE_SET.TIME_SERIES
}
interface HeatMapItem extends DisplayDataItem {
  dataType: typeof DATA_TYPE_SET.HEAT_MAP
}
