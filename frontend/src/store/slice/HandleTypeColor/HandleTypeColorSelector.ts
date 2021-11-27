import { RootState } from '../../store'

export const handleTypeColorSelector = (key: string) => (state: RootState) => {
  if (state.handleColor.colorMap[key] != null) {
    return state.handleColor.colorMap[key]
  } else {
    return undefined
  }
}

export const nextColorKeySelector = (state: RootState) =>
  state.handleColor.nextKey
