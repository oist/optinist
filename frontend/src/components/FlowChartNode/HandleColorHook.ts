import React from 'react'
import { useDispatch, useSelector } from 'react-redux'
import { selectHandleTypeColor } from 'store/slice/HandleTypeColor/HandleTypeColorSelectors'
import { addColor } from 'store/slice/HandleTypeColor/HandleTypeColorSlice'

export function useHandleColor(type: string) {
  const dispatch = useDispatch()
  const color = useSelector(selectHandleTypeColor(type))
  React.useEffect(() => {
    if (color === undefined) {
      dispatch(addColor(type))
    }
  }, [type, color, dispatch])
  return color
}
