export const PLOT_DATA_SLICE_NAME = 'plotData'

/**
 * ${nodeId}/${outputKey}
 */
export type PlotDataKey = `${string}/${string}`

export type PlotData = {
  timeSeriesDataMap: {
    // keyの方はPlotDataKeyを使いたいけど、テンプレートリテラルはオブジェクトのkeyとして使えない...
    [key: string]: TimeSeriesData
  }
  heatMapDataMap: {
    [key: string]: HeatMapData
  }
}

export type TimeSeriesData = {
  data: {
    [key: string]: {
      [key: number]: number
    }
  }
}

export type HeatMapData = {
  data: number[][]
}
