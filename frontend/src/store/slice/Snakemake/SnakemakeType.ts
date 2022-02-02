export const Snakemake_SLICE_NAME = 'Snakemake'

export type SnakemakeType = {
  SnakemakeList: SnakemakeListType
}

export type SnakemakeListDTO = {
  [name: string]: { children: SnakemakeListDTO }
}

export type SnakemakeListType = {
  [SnakemakeName: string]: SnakemakeNodeType
}

export type SnakemakeNodeType = SnakemakeChild | SnakemakeParent

export type SnakemakeChild = {
  type: 'child'
  value: unknown
  path: string
}

export type SnakemakeParent = {
  type: 'parent'
  children: {
    [name: string]: SnakemakeNodeType
  }
}
