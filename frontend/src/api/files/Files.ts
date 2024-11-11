import { AxiosProgressEvent } from "axios"

import { BASE_URL } from "const/API"
import axios from "utils/axios"

export const FILE_TREE_TYPE_SET = {
  IMAGE: "image",
  CSV: "csv",
  HDF5: "hdf5",
  FLUO: "fluo",
  BEHAVIOR: "behavior",
  MATLAB: "matlab",
  MICROSCOPE: "microscope",
  ALL: "all",
} as const

export type FILE_TREE_TYPE =
  (typeof FILE_TREE_TYPE_SET)[keyof typeof FILE_TREE_TYPE_SET]

export type TreeNodeTypeDTO = DirNodeDTO | FileNodeDTO

export interface NodeBaseDTO {
  path: string
  name: string
  isdir: boolean
  shape: []
}

export interface DirNodeDTO extends NodeBaseDTO {
  isdir: true
  nodes: TreeNodeTypeDTO[]
}

export interface FileNodeDTO extends NodeBaseDTO {
  isdir: false
}

export type GetStatusViaUrl = {
  total: number
  current: number
  error: string | null
}

export async function getFilesTreeApi(
  workspaceId: number,
  fileType: FILE_TREE_TYPE,
): Promise<TreeNodeTypeDTO[]> {
  const response = await axios.get(`${BASE_URL}/files/${workspaceId}`, {
    params: {
      file_type: fileType,
    },
  })
  return response.data
}

export async function uploadFileApi(
  workspaceId: number,
  fileName: string,
  config: {
    onUploadProgress: (progressEvent: AxiosProgressEvent) => void
  },
  formData: FormData,
): Promise<{ file_path: string }> {
  const response = await axios.post(
    `${BASE_URL}/files/${workspaceId}/upload/${fileName}`,
    formData,
    config,
  )
  return response.data
}

export async function deleteFileTreeApi(
  workspaceId: number,
  fileName: string,
): Promise<boolean> {
  const response = await axios.delete(
    `${BASE_URL}/${workspaceId}/delete/${fileName}`,
  )
  return response.data
}

export async function updateShapeApi(
  workspaceId: number,
  fileName: string,
): Promise<boolean> {
  const response = await axios.post(
    `${BASE_URL}/files/${workspaceId}/shape/${fileName}`,
  )
  return response.data
}

export const uploadViaUrlApi = async (
  workspaceId: number,
  url: string,
): Promise<{ file_name: string }> => {
  const res = await axios.post(
    `${BASE_URL}/files/${workspaceId}/download`,
    { url },
    // config,
  )
  return res.data
}

export const getStatusLoadViaUrlApi = async (
  workspaceId: number,
  file_name: string,
): Promise<GetStatusViaUrl> => {
  const res = await axios.get(
    `${BASE_URL}/files/${workspaceId}/download/status?file_name=${file_name}`,
  )
  return res.data
}
