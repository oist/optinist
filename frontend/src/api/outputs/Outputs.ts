import axios from "utils/axios"
import { PlotMetaData } from "store/slice/DisplayData/DisplayDataType"
import { BASE_URL } from "const/API"

export type TimeSeriesData = {
  [key: string]: {
    [key: string]: number
  }
}

export async function getTimeSeriesInitDataApi(path: string): Promise<{
  data: TimeSeriesData
  xrange: string[]
  std: TimeSeriesData
  meta?: PlotMetaData
}> {
  const response = await axios.get(`${BASE_URL}/outputs/inittimedata/${path}`)
  return response.data
}

export async function getTimeSeriesDataByIdApi(
  path: string,
  index: string,
): Promise<{
  data: TimeSeriesData
  xrange: string[]
  std: TimeSeriesData
  meta?: PlotMetaData
}> {
  const response = await axios.get(`${BASE_URL}/outputs/timedata/${path}`, {
    params: {
      index: index,
    },
  })
  return response.data
}

export async function getTimeSeriesAllDataApi(path: string): Promise<{
  data: TimeSeriesData
  xrange: string[]
  std: TimeSeriesData
  meta?: PlotMetaData
}> {
  const response = await axios.get(`${BASE_URL}/outputs/alltimedata/${path}`)
  return response.data
}

export type HeatMapData = number[][]

export async function getHeatMapDataApi(path: string): Promise<{
  data: HeatMapData
  columns: string[]
  index: string[]
  meta?: PlotMetaData
}> {
  const response = await axios.get(`${BASE_URL}/outputs/data/${path}`)
  return response.data
}

export type ImageData = number[][][]

export async function getImageDataApi(
  path: string,
  params: {
    workspaceId: number
    startIndex?: number
    endIndex?: number
  },
): Promise<{ data: ImageData; meta?: PlotMetaData }> {
  const response = await axios.get(`${BASE_URL}/outputs/image/${path}`, {
    params: {
      workspace_id: params.workspaceId,
      start_index: params.startIndex,
      end_index: params.endIndex,
    },
  })
  return response.data
}

export type CsvData = number[][]

export async function getCsvDataApi(
  path: string,
  params: { workspaceId: number },
): Promise<{ data: CsvData; meta?: PlotMetaData }> {
  const response = await axios.get(`${BASE_URL}/outputs/csv/${path}`, {
    params: { workspace_id: params.workspaceId },
  })

  return response.data
}

export type RoiData = number[][][]

export async function getRoiDataApi(
  path: string,
  params: { workspaceId: number },
): Promise<{ data: RoiData; meta?: PlotMetaData }> {
  const response = await axios.get(`${BASE_URL}/outputs/image/${path}`, {
    params: { workspace_id: params.workspaceId },
  })
  return response.data
}

export type ScatterData = {
  [key: string]: {
    [key: number]: number
  }
}

export async function getScatterDataApi(
  path: string,
): Promise<{ data: ScatterData; meta?: PlotMetaData }> {
  const response = await axios.get(`${BASE_URL}/outputs/data/${path}`, {})
  return response.data
}

export type BarData = {
  [key: string]: {
    [key: number]: number
  }
}

export async function getBarDataApi(path: string): Promise<{
  data: BarData
  columns: string[]
  index: string[]
  meta?: PlotMetaData
}> {
  const response = await axios.get(`${BASE_URL}/outputs/data/${path}`, {})
  return response.data
}

export type HTMLData = string

export async function getHTMLDataApi(
  path: string,
): Promise<{ data: HTMLData; meta?: PlotMetaData }> {
  const response = await axios.get(`${BASE_URL}/outputs/html/${path}`, {})
  return response.data
}

export async function addRoiApi(
  path: string,
  data: { posx: number; posy: number; sizex: number; sizey: number },
): Promise<{ data: HTMLData; meta?: PlotMetaData }> {
  const response = await axios.post(
    `${BASE_URL}/outputs/image/${path}/add_roi`,
    data,
  )
  return response.data
}

export async function mergeRoiApi(
  path: string,
  data: { ids: number[] },
): Promise<{ data: HTMLData; meta?: PlotMetaData }> {
  const response = await axios.post(
    `${BASE_URL}/outputs/image/${path}/merge_roi`,
    data,
  )
  return response.data
}

export async function deleteRoiApi(
  path: string,
  data: { ids: number[] },
): Promise<{ data: HTMLData; meta?: PlotMetaData }> {
  const response = await axios.post(
    `${BASE_URL}/outputs/image/${path}/delete_roi`,
    data,
  )
  return response.data
}

export type LineData = number[][]

export async function getLineDataApi(
  path: string,
): Promise<{ data: LineData; columns: number[]; index: number[] }> {
  const response = await axios.get(`${BASE_URL}/outputs/data/${path}`, {})
  return response.data
}

export type PolarData = number[][]

export async function getPolarDataApi(
  path: string,
): Promise<{ data: PolarData; columns: number[]; index: number[] }> {
  const response = await axios.get(`${BASE_URL}/outputs/data/${path}`, {})
  return response.data
}

export type HistogramData = number[][]

export async function getHistogramDataApi(
  path: string,
): Promise<{ data: HistogramData }> {
  const response = await axios.get(`${BASE_URL}/outputs/data/${path}`, {})
  return response.data
}

export type PieData = number[][]

export async function getPieDataApi(
  path: string,
): Promise<{ data: PieData; columns: string[] }> {
  const response = await axios.get(`${BASE_URL}/outputs/data/${path}`, {})
  return response.data
}
