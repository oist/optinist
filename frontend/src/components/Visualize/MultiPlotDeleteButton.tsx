import React from 'react'
import { useDispatch } from 'react-redux'

import IconButton from '@mui/material/IconButton'
import CloseIcon from '@mui/icons-material/Close'
import { deleteVisualizeItem } from 'store/slice/VisualizeItem/VisualizeItemSlice'

export const MultiPlotDeleteButton = React.memo<{ itemId: number }>(
  ({ itemId }) => {
    const dispatch = useDispatch()

    const onClick: React.MouseEventHandler<HTMLButtonElement> = (e) => {
      e.stopPropagation() // 親divのonClickを反応させないため
      dispatch(deleteVisualizeItem(itemId))
    }
    return (
      <IconButton onClick={onClick} size="large">
        <CloseIcon />
      </IconButton>
    )
  },
)
