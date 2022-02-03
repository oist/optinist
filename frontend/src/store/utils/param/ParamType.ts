export type ParamMap = {
  [paramKey: string]: ParamType
}

export type ParamType = ParamParent | ParamChild

export type ParamParent = {
  type: 'parent'
  children: {
    [key: string]: ParamType
  }
}

export type ParamChild = {
  type: 'child'
  value: unknown
  path: string
}

export type ParamDTO = {
  [key: string]: unknown
}
