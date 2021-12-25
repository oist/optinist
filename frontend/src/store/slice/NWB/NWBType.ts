export const NWB_SLICE_NAME = 'NWB'

export type NWBType = {
  nwbList: NWBListType
}

export type NWBListDTO = {
  [name: string]: { children: NWBListDTO }
}

export type NWBListType = {
  [nwbName: string]: NWBNodeType
}

export type NWBNodeType = NWBChild | NWBParent

export type NWBChild = {
  type: 'child'
  value: unknown
  path: string
}

export type NWBParent = {
  type: 'parent'
  children: {
    [name: string]: NWBNodeType
  }
}
