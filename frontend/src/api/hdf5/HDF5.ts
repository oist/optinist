import axios from 'utils/axios'

import { BASE_URL } from 'const/API'

export type HDF5TreeDTO = HDF5DirDTO | HDF5FileDTO

export interface HDF5DirDTO {
  isDir: true
  name: string
  nodes: HDF5TreeDTO[]
  path: string
}

export interface HDF5FileDTO {
  isDir: false
  name: string
  shape: [number]
  path: string
  nbytes: string
}

export async function getHDF5TreeApi(path: string): Promise<HDF5TreeDTO[]> {
  const response = await axios.get(`${BASE_URL}/hdf5/${path}`)
  return response.data
}
