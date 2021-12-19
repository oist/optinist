export const NWB_SLICE_NAME = 'NWB'

export type NWBType = {
  params: ParamType
  nwbList: NWBListType
}

export type ParamType = {
  [name: string]: unknown
}

export type NWBListDTO = {
  [name: string]: { children: NWBListDTO }
}

export type NWBListType = {
  [nwbName: string]: NWBNodeType
}

export type NWBNodeType = NWBChild | NWBParent

type NWBChild = {
  type: 'child'
  value: unknown
}

type NWBParent = {
  type: 'parent'
  children: {
    [name: string]: NWBNodeType
  }
}
