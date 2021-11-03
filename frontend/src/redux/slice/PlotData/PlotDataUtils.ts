import { PlotDataKey } from './PlotDataType'

export function toPlotDataKey(nodeId: string, outputKey: string): PlotDataKey {
  return `${nodeId}/${outputKey}`
}
