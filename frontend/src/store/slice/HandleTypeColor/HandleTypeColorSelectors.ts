import { RootState } from '../../store'

export const selectHandleTypeColor = (key: string) => (state: RootState) => {
  if (state.handleColor.colorMap[key] != null) {
    return state.handleColor.colorMap[key]
  } else {
    return undefined
  }
}
