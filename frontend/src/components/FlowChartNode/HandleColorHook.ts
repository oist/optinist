import React from 'react'
import { useDispatch, useSelector } from 'react-redux'
import { handleTypeColorSelector } from 'store/slice/HandleTypeColor/HandleTypeColorSelector'
import { addColor } from 'store/slice/HandleTypeColor/HandleTypeColor'

export function useHandleColor(type: string) {
  const dispatch = useDispatch()
  const color = useSelector(handleTypeColorSelector(type))
  React.useEffect(() => {
    if (color === undefined) {
      dispatch(addColor(type))
    }
  }, [type, color, dispatch])
  return color
}
