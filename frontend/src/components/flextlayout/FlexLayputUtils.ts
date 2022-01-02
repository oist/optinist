import { TAB_COMPONENT_TYPE } from 'store/slice/LayoutTab/LayoutTabType'

export function toLayoutTab(
  tabId: string,
  tabComponentType: TAB_COMPONENT_TYPE,
  name: string,
) {
  return {
    id: tabId,
    component: tabComponentType,
    name,
    type: 'tab',
  }
}
