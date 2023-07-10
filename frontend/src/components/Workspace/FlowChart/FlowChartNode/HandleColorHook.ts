import { useSelector } from 'react-redux'
import { selectHandleTypeColor } from 'store/slice/HandleTypeColor/HandleTypeColorSelectors'

export function useHandleColor(type: string) {
  const color = useSelector(selectHandleTypeColor(type))
  return color
}
