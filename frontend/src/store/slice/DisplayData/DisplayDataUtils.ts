import { DATA_TYPE, DATA_TYPE_SET } from './DisplayDataType'

/**
 * サーバーサイドの値をマッピング
 */
export function toDataType(value: string): DATA_TYPE {
  switch (value) {
    case 'images':
      return DATA_TYPE_SET.IMAGE
    case 'timeseries':
      return DATA_TYPE_SET.TIME_SERIES
    case 'heatmap':
      return DATA_TYPE_SET.HEAT_MAP
    case 'roi':
      return DATA_TYPE_SET.ROI
    case 'scatter':
      return DATA_TYPE_SET.SCATTER
    case 'bar':
      return DATA_TYPE_SET.BAR
    case 'html':
      return DATA_TYPE_SET.HTML
    case 'histogram':
      return DATA_TYPE_SET.HISTOGRAM
    case 'line':
      return DATA_TYPE_SET.LINE
    case 'pie':
      return DATA_TYPE_SET.PIE
    case 'polar':
      return DATA_TYPE_SET.POLAR
    default:
      throw new Error(`failed to map dataType: ${value}`)
  }
}
