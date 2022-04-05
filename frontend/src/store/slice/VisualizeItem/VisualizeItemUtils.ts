import { DATA_TYPE_SET } from '../DisplayData/DisplayDataType'
import {
  VisualaizeItemType,
  DisplayDataItem,
  VISUALIZE_ITEM_TYPE_SET,
  ImageItem,
  TimeSeriesItem,
  HeatMapItem,
  RoiItem,
  ScatterItem,
  CsvItem,
} from './VisualizeItemType'

export function isDisplayDataItem(
  item: VisualaizeItemType,
): item is DisplayDataItem {
  return item.itemType === VISUALIZE_ITEM_TYPE_SET.DISPLAY_DATA
}

export function isImageItem(item: VisualaizeItemType): item is ImageItem {
  return (
    item.itemType === VISUALIZE_ITEM_TYPE_SET.DISPLAY_DATA &&
    item.dataType === DATA_TYPE_SET.IMAGE
  )
}

export function isTimeSeriesItem(
  item: VisualaizeItemType,
): item is TimeSeriesItem {
  return (
    item.itemType === VISUALIZE_ITEM_TYPE_SET.DISPLAY_DATA &&
    item.dataType === DATA_TYPE_SET.TIME_SERIES
  )
}

export function isCsvItem(item: VisualaizeItemType): item is CsvItem {
  return (
    item.itemType === VISUALIZE_ITEM_TYPE_SET.DISPLAY_DATA &&
    item.dataType === DATA_TYPE_SET.CSV
  )
}

export function isHeatMapItem(item: VisualaizeItemType): item is HeatMapItem {
  return (
    item.itemType === VISUALIZE_ITEM_TYPE_SET.DISPLAY_DATA &&
    item.dataType === DATA_TYPE_SET.HEAT_MAP
  )
}

export function isRoiItem(item: VisualaizeItemType): item is RoiItem {
  return (
    item.itemType === VISUALIZE_ITEM_TYPE_SET.DISPLAY_DATA &&
    item.dataType === DATA_TYPE_SET.ROI
  )
}

export function isScatterItem(item: VisualaizeItemType): item is ScatterItem {
  return (
    item.itemType === VISUALIZE_ITEM_TYPE_SET.DISPLAY_DATA &&
    item.dataType === DATA_TYPE_SET.SCATTER
  )
}
