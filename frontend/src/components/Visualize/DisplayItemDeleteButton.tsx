import React from 'react'
import { useDispatch, useSelector } from 'react-redux'

import IconButton from '@material-ui/core/IconButton'
import CloseIcon from '@material-ui/icons/Close'
import { deleteVisualizeItem } from 'store/slice/VisualizeItem/VisualizeItemSlice'
import {
  selectVisualizeDataFilePath,
  selectVisualizeDataType,
} from 'store/slice/VisualizeItem/VisualizeItemSelectors'
import { deleteDisplayItem } from 'store/slice/DisplayData/DisplayDataSlice'

export const DisplayItemDeleteButton = React.memo<{ itemId: number }>(
  ({ itemId }) => {
    const dispatch = useDispatch()

    const dataType = useSelector(selectVisualizeDataType(itemId))
    const filePath = useSelector(selectVisualizeDataFilePath(itemId))

    const onClick: React.MouseEventHandler<HTMLButtonElement> = (e) => {
      e.stopPropagation() // 親divのonClickを反応させないため
      dispatch(deleteVisualizeItem(itemId))
      // visualize Itemで同じpathのデータ個数を調べて、1だったら、displayも削除
      dispatch(deleteDisplayItem({ dataType, filePath }))
    }
    return (
      <IconButton onClick={onClick}>
        <CloseIcon />
      </IconButton>
    )
  },
)
