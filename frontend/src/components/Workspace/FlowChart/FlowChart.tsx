import React from 'react'
import { useSelector } from 'react-redux'
import { styled } from '@mui/material/styles'
import { DndProvider } from 'react-dnd'
import { HTML5Backend } from 'react-dnd-html5-backend'

import { AlgorithmTreeView } from './TreeView'
import { ReactFlowComponent } from './ReactFlowComponent'
import RightDrawer from './RightDrawer'
import { selectRightDrawerIsOpen } from 'store/slice/RightDrawer/RightDrawerSelectors'
import { UseRunPipelineReturnType } from 'store/slice/Pipeline/PipelineHook'
import { CurrentPipelineInfo } from 'components/common/CurrentPipelineInfo'
import { DRAWER_WIDTH, RIGHT_DRAWER_WIDTH } from 'const/Layout'
import { Box } from '@mui/material'
import { grey } from '@mui/material/colors'

const FlowChart = React.memo<UseRunPipelineReturnType>((props) => {
  const open = useSelector(selectRightDrawerIsOpen)
  return (
    <RootDiv>
      <DndProvider backend={HTML5Backend}>
        <Box
          sx={{
            width: DRAWER_WIDTH,
          }}
          borderRight={1}
          borderColor={grey[300]}
        >
          <CurrentPipelineInfo />
          <DrawerContents>
            <AlgorithmTreeView />
          </DrawerContents>
        </Box>

        <MainContents open={open}>
          <ReactFlowComponent {...props} />
        </MainContents>
      </DndProvider>
      <RightDrawer />
    </RootDiv>
  )
})

const RootDiv = styled('div')({
  display: 'flex',
})

const DrawerContents = styled('div')({
  overflow: 'auto',
})

const MainContents = styled('main')<{ open: boolean }>(
  ({ theme }) => ({
    flexDirection: 'column',
    flexGrow: 1,
    height: '90vh',
    transition: theme.transitions.create('margin', {
      easing: theme.transitions.easing.sharp,
      duration: theme.transitions.duration.leavingScreen,
    }),
    marginRight: -RIGHT_DRAWER_WIDTH,
  }),
  ({ open, theme }) =>
    open
      ? {
          transition: theme.transitions.create('margin', {
            easing: theme.transitions.easing.easeOut,
            duration: theme.transitions.duration.enteringScreen,
          }),
          marginRight: 0,
        }
      : undefined,
)

export default FlowChart
