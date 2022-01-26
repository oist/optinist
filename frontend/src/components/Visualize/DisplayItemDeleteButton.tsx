import React from 'react'
import { useDispatch, useSelector } from 'react-redux'

import IconButton from '@material-ui/core/IconButton'
import CloseIcon from '@material-ui/icons/Close'
import { deleteVisualizeItem } from 'store/slice/VisualizeItem/VisualizeItemSlice'
import {
  selectVisualizeDataFilePath,
  selectVisualizeDataType,
} from 'store/slice/VisualizeItem/VisualizeItemSelectors'
import {
  deleteDisplayHeatMapItem,
  deleteDisplayImageItem,
  deleteDisplayRoiItem,
  deleteDisplayScatterItem,
  deleteDisplayTableItem,
  deleteDisplayTimeSeriesItem,
} from 'store/slice/DisplayData/DisplayDataSlice'

export const DisplayItemDeleteButton = React.memo<{ itemId: number }>(
  ({ itemId }) => {
    const dispatch = useDispatch()

    const dataType = useSelector(selectVisualizeDataType(itemId))
    const filePath = useSelector(selectVisualizeDataFilePath(itemId))

    const onClick: React.MouseEventHandler<HTMLButtonElement> = (e) => {
      e.stopPropagation() // 親divのonClickを反応させないため
      dispatch(deleteVisualizeItem(itemId))
      // visualize Itemで同じpathのデータ個数を調べて、1だったら、displayも削除
      if (filePath !== null) {
        if (dataType === 'image') {
          dispatch(deleteDisplayImageItem({ filePath }))
        } else if (dataType === 'timeSeries') {
          dispatch(deleteDisplayTimeSeriesItem({ filePath }))
        } else if (dataType === 'table') {
          dispatch(deleteDisplayTableItem({ filePath }))
        } else if (dataType === 'heatMap') {
          dispatch(deleteDisplayHeatMapItem({ filePath }))
        } else if (dataType === 'roi') {
          dispatch(deleteDisplayRoiItem({ filePath }))
        } else if (dataType === 'scatter') {
          dispatch(deleteDisplayScatterItem({ filePath }))
        } else {
          throw new Error('invalid item Type')
        }
      }
    }
    return (
      <IconButton onClick={onClick}>
        <CloseIcon />
      </IconButton>
    )
  },
)
