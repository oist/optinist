import { BASE_URL } from "const/API"
import axios from "utils/axios"

export type MatlabTreeDTO = MatlabDirDTO | MatlabFileDTO

export interface MatlabDirDTO {
  isDir: true
  name: string
  nodes: MatlabTreeDTO[]
  path: string
  dataType: string | null
}

export interface MatlabFileDTO {
  isDir: false
  name: string
  shape: [number]
  path: string
  nbytes: string
  dataType: string | null
}

export async function getMatlabTreeApi(
  path: string,
  workspaceId: number,
): Promise<MatlabTreeDTO[]> {
  const response = await axios.get(
    `${BASE_URL}/mat/${path}?workspace_id=${workspaceId}`,
  )
  return response.data
}
