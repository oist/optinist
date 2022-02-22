export const HDF5_SLICE_NAME = 'hdf5'

export type HDF5TreeDTO = HDF5DirDTO | HDF5FileDTO

export interface HDF5DirDTO {
  isDir: true
  name: string
  nodes: HDF5TreeDTO[]
}

export interface HDF5FileDTO {
  isDir: false
  name: string
  shape: [number]
  path: string
}

export interface HDF5Tree {
  isLatest: boolean
  isLoading: boolean
  tree: HDF5TreeDTO[]
}

export const HDF5_TYPE_SET = {
  ALL: 'all',
} as const

export type HDF5_TYPE = typeof HDF5_TYPE_SET[keyof typeof HDF5_TYPE_SET]
