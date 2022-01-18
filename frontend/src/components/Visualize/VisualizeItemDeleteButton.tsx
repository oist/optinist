import React from 'react'
import { useDispatch } from 'react-redux'

import IconButton from '@material-ui/core/IconButton'
import CloseIcon from '@material-ui/icons/Close'
import { deleteItem } from 'store/slice/VisualizeItem/VisualizeItemSlice'

export const VisualizeItemDeleteButton = React.memo<{ itemId: number }>(
  ({ itemId }) => {
    const dispatch = useDispatch()
    const onClick: React.MouseEventHandler<HTMLButtonElement> = (e) => {
      e.stopPropagation() // 親divのonClickを反応させないため
      dispatch(deleteItem(itemId))
    }
    return (
      <IconButton onClick={onClick}>
        <CloseIcon />
      </IconButton>
    )
  },
)
