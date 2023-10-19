import axios from "utils/axios"
import qs from "qs"
import {
  ItemsWorkspace,
  WorkspaceDataDTO,
  ListUserShareWorkspaceDTO,
} from "store/slice/Workspace/WorkspaceType"

export type WorkspacePostDataDTO = { name: string; id?: number }

export const getWorkspaceApi = async (id: number): Promise<ItemsWorkspace> => {
  const response = await axios.get(`/workspace/${id}`)
  return response.data
}

export const getWorkspacesApi = async (params: {
  [key: string]: number
}): Promise<WorkspaceDataDTO> => {
  const paramsNew = qs.stringify(params, { indices: false })
  const response = await axios.get(`/workspaces?${paramsNew}`)
  return response.data
}

export const delWorkspaceApi = async (id: number): Promise<boolean> => {
  const response = await axios.delete(`/workspace/${id}`)
  return response.data
}

export const postWorkspaceApi = async (
  data: WorkspacePostDataDTO,
): Promise<ItemsWorkspace> => {
  const response = await axios.post(`/workspace`, data)
  return response.data
}

export const putWorkspaceApi = async (
  data: WorkspacePostDataDTO,
): Promise<ItemsWorkspace> => {
  const response = await axios.put(`/workspace/${data.id}`, { name: data.name })
  return response.data
}

export const getListUserShareWorkspaceApi = async (
  id: number,
): Promise<ListUserShareWorkspaceDTO> => {
  const response = await axios.get(`/workspace/share/${id}/status`)
  return response.data
}

export const postListUserShareWorkspaceApi = async (
  id: number,
  data: { user_ids: number[] },
): Promise<boolean> => {
  const response = await axios.post(`/workspace/share/${id}/status`, data)
  return response.data
}
