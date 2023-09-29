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
  BarItem,
  HistogramItem,
  LineItem,
  PieItem,
  PolarItem,
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

export function isBarItem(item: VisualaizeItemType): item is BarItem {
  return (
    item.itemType === VISUALIZE_ITEM_TYPE_SET.DISPLAY_DATA &&
    item.dataType === DATA_TYPE_SET.BAR
  )
}

export function isHistogramItem(
  item: VisualaizeItemType,
): item is HistogramItem {
  return (
    item.itemType === VISUALIZE_ITEM_TYPE_SET.DISPLAY_DATA &&
    item.dataType === DATA_TYPE_SET.HISTOGRAM
  )
}

export function isLineItem(item: VisualaizeItemType): item is LineItem {
  return (
    item.itemType === VISUALIZE_ITEM_TYPE_SET.DISPLAY_DATA &&
    item.dataType === DATA_TYPE_SET.LINE
  )
}

export function isPieItem(item: VisualaizeItemType): item is PieItem {
  return (
    item.itemType === VISUALIZE_ITEM_TYPE_SET.DISPLAY_DATA &&
    item.dataType === DATA_TYPE_SET.PIE
  )
}

export function isPolarItem(item: VisualaizeItemType): item is PolarItem {
  return (
    item.itemType === VISUALIZE_ITEM_TYPE_SET.DISPLAY_DATA &&
    item.dataType === DATA_TYPE_SET.POLAR
  )
}
