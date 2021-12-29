import { DATA_TYPE } from '../DisplayData/DisplayDataType'
import {
  TabType,
  TAB_COMPONENT_TYPE_SET,
  DisplayTab,
  ParamFormTab,
} from './LayoutTabType'

export function isDisplayTab(tab: TabType): tab is DisplayTab {
  return tab.type === TAB_COMPONENT_TYPE_SET.DISPLAY_DATA
}

export function isParamFormTab(tab: TabType): tab is ParamFormTab {
  return tab.type === TAB_COMPONENT_TYPE_SET.PARAM_FORM
}

export function equalsDisplayTabInfo(
  a: TabType,
  b: { nodeId: string; dataType: DATA_TYPE; filePath: string },
) {
  if (isDisplayTab(a)) {
    return (
      a.nodeId === b.nodeId &&
      a.dataType === b.dataType &&
      a.filePath === b.filePath
    )
  } else {
    return false
  }
}
