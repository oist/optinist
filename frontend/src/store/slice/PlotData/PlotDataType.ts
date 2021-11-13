export const PLOT_DATA_SLICE_NAME = 'plotData'

/**
 * ${nodeId}/${outputKey}
 */
export type PlotDataKey = `${string}/${string}`

export type PlotData = {
  timeSeriesDataMap: {
    // keyの方はPlotDataKeyを使いたいけど、テンプレートリテラルはオブジェクトのkeyとして使えない...
    [path: string]: {
      data: TimeSeriesData
    }
  }
  heatMapDataMap: {
    [path: string]: {
      data: HeatMapData
    }
  }
  imageDataMap: {
    // PlotDataKeyの時と、ファイル選択ノードのnodeIdの時がある
    [path: string]: {
      activeIndex: number
      data: ImageData
    }
  }
}

export type TimeSeriesData = {
  [key: string]: {
    [key: number]: number
  }
}

export type HeatMapData = number[][]

export type ImageData = number[][][]
