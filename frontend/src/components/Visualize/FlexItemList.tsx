import React from 'react'
import { useSelector } from 'react-redux'

import Box from '@mui/material/Box'

import { twoDimarrayEqualityFn } from 'utils/EqualityUtils'

import { selectVisualizeItemLayout } from 'store/slice/VisualizeItem/VisualizeItemSelectors'

import { VisualizeItemAddButton } from './VisualizeItemAddButton'
import { VisualizeItem } from './VisualizeItem'

export const FlexItemList: React.FC = () => {
  const layout = useSelector(selectVisualizeItemLayout, twoDimarrayEqualityFn)
  return (
    <Box display="flex" flexWrap="wrap" flexDirection="column" p={1} m={1}>
      {layout.map((row) => (
        <Box display="flex" flexDirection="row">
          {row.map((itemId) => (
            <VisualizeItem itemId={itemId} key={itemId} />
          ))}
        </Box>
      ))}
      <VisualizeItemAddButton />
    </Box>
  )
}
