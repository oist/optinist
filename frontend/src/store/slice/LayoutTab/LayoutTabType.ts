import { DATA_TYPE } from '../DisplayData/DisplayDataType'

export const LAYOUT_TAB_SLICE_NAME = 'layoutTab'

export type LayoutTab = {
  tabs: {
    [tabId: string]: TabType
  }
}

export type TabType = ParamFormTab | DisplayTab

export const TAB_COMPONENT_TYPE_SET = {
  DISPLAY_DATA: 'displayData',
  PARAM_FORM: 'paramForm',
} as const

export type TAB_COMPONENT_TYPE =
  typeof TAB_COMPONENT_TYPE_SET[keyof typeof TAB_COMPONENT_TYPE_SET]

interface TabBaseType<T extends TAB_COMPONENT_TYPE> {
  tabId: string
  type: T
  nodeId: string
}

export interface DisplayTab extends TabBaseType<'displayData'> {
  filePath: string
  dataType: DATA_TYPE
}

export interface ParamFormTab extends TabBaseType<'paramForm'> {}
