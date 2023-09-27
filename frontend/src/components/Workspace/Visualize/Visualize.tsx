import React from 'react'
import { styled } from '@mui/material/styles'
import { FlexItemList } from './FlexItemList'
import { VisualizeItemEditor } from './VisualizeItemEditor'
import { CurrentPipelineInfo } from 'components/common/CurrentPipelineInfo'
import { CONTENT_HEIGHT, DRAWER_WIDTH } from 'const/Layout'
import { Box } from '@mui/material'
import { grey } from '@mui/material/colors'

const Visualize: React.FC = () => {
  return (
    <RootDiv>
      <Box
        sx={{
          width: DRAWER_WIDTH,
        }}
        borderRight={1}
        borderColor={grey[300]}
      >
        <CurrentPipelineInfo />
        <StyledDrawerContents>
          <VisualizeItemEditor />
        </StyledDrawerContents>
      </Box>

      <MainContents>
        <FlexItemList />
      </MainContents>
    </RootDiv>
  )
}

const RootDiv = styled('div')({
  display: 'flex',
})

const StyledDrawerContents = styled('div')({
  overflow: 'auto',
})

const MainContents = styled('main')({
  display: 'flex',
  flexDirection: 'column',
  flexGrow: 1,
  minHeight: CONTENT_HEIGHT
})

export default Visualize
