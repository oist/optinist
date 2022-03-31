import React, { useEffect, useState } from 'react'
import { useDispatch, useSelector } from 'react-redux'

import { useTheme } from '@mui/material/styles'
import Box from '@mui/material/Box'
import Paper from '@mui/material/Paper'

import { arrayEqualityFn } from 'utils/EqualityUtils'

import {
  selectSelectedVisualizeItemId,
  selectVisualizeDataFilePath,
  selectVisualizeDataNodeId,
  selectVisualizeDataType,
  selectVisualizeItemIdList,
  selectVisualizeItemType,
} from 'store/slice/VisualizeItem/VisualizeItemSelectors'

import { VisualizeItemAddButton } from './VisualizeItemAddButton'
import { DisplayItemDeleteButton } from './DisplayItemDeleteButton'
import { MultiPlotItem } from './MultiPlotItem'
import { DisplayDataItem } from './DisplayDataItem'
import {
  selectItem,
  setDisplayDataPath,
} from 'store/slice/VisualizeItem/VisualizeItemSlice'
import { RootState } from 'store/store'
import { MultiPlotDeleteButton } from './MultiPlotDeleteButton'
import { deleteDisplayItem } from 'store/slice/DisplayData/DisplayDataSlice'
import { DATA_TYPE } from 'store/slice/DisplayData/DisplayDataType'
import { FilePathSelect } from './FilePathSelect'
import { Grid } from '@mui/material'

export const VisualizeItems: React.FC = () => {
  return <FlexItemList />
}

const FlexItemList: React.FC = () => {
  const itemIdList = useSelector(selectVisualizeItemIdList, arrayEqualityFn)
  return (
    <Box display="flex" flexWrap="wrap" p={1} m={1}>
      {itemIdList.map((itemId) => (
        <RowItem itemId={itemId} key={itemId} />
      ))}
      <VisualizeItemAddButton />
    </Box>
  )
}

const RowItem = React.memo<{ itemId: number }>(({ itemId }) => {
  const dispatch = useDispatch()
  const onClick = () => {
    dispatch(selectItem(itemId))
  }
  const isSelected = useSelector(
    (state: RootState) => selectSelectedVisualizeItemId(state) === itemId,
  )
  const theme = useTheme()

  return (
    <Paper
      variant="outlined"
      key={itemId}
      style={{
        width: '100%',
        margin: theme.spacing(1),
        padding: theme.spacing(1),
        cursor: 'pointer',
        borderColor: isSelected ? theme.palette.primary.light : undefined,
      }}
      onClick={onClick}
    >
      <Item itemId={itemId} />
    </Paper>
  )
})

const Item = React.memo<{ itemId: number }>(({ itemId }) => {
  const dispatch = useDispatch()
  const itemType = useSelector(selectVisualizeItemType(itemId))
  const filePath = useSelector(selectVisualizeDataFilePath(itemId))
  const nodeId = useSelector(selectVisualizeDataNodeId(itemId))
  const dataType = useSelector(selectVisualizeDataType(itemId))

  const [prevItem, setPrevItem] = useState<{
    dataType: DATA_TYPE
    filePath: string | null
  }>({
    dataType: 'image',
    filePath: null,
  })

  useEffect(() => {
    setPrevItem({ dataType, filePath })
  }, [filePath, dataType])

  const onSelect = (nodeId: string, filePath: string, dataType: DATA_TYPE) => {
    dispatch(setDisplayDataPath({ itemId, nodeId, filePath, dataType }))
    dispatch(deleteDisplayItem(prevItem))
  }

  return (
    <Box>
      <Grid container spacing={2}>
        <Grid item xs={10}>
          <FilePathSelect
            selectedNodeId={nodeId}
            selectedFilePath={filePath}
            onSelect={onSelect}
          />
        </Grid>
        <Grid item xs={2} display="flex" justifyContent="flex-end">
          <DeleteButton itemType={itemType} itemId={itemId} />
        </Grid>
      </Grid>
      <ItemByType itemType={itemType} itemId={itemId} />
    </Box>
  )
})

const DeleteButton: React.FC<{
  itemType: string
  itemId: number
}> = ({ itemType, itemId }) => {
  switch (itemType) {
    case 'displayData':
      return <DisplayItemDeleteButton itemId={itemId} />
    case 'MultiPlot':
      return <MultiPlotDeleteButton itemId={itemId} />
    default:
      throw new Error('itemType Error')
  }
}

const ItemByType = React.memo<{
  itemType: string
  itemId: number
}>(({ itemType, itemId }) => {
  switch (itemType) {
    case 'MultiPlot':
      return <MultiPlotItem itemId={itemId} />
    case 'displayData':
      return <DisplayDataItem itemId={itemId} />
    default:
      throw new Error('itemType Error')
  }
})
