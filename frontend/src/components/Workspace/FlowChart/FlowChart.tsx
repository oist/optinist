import React from 'react'
import { useSelector } from 'react-redux'
import Drawer, { drawerClasses } from '@mui/material/Drawer'
import { default as MuiToolbar } from '@mui/material/Toolbar'
import { styled } from '@mui/material/styles'
import { DndProvider } from 'react-dnd'
import { HTML5Backend } from 'react-dnd-html5-backend'

import { AlgorithmTreeView } from './TreeView'
import { ReactFlowComponent } from './ReactFlowComponent'
import RightDrawer, { rightDrawerWidth } from './RightDrawer'
import { selectRightDrawerIsOpen } from 'store/slice/RightDrawer/RightDrawerSelectors'
import { UseRunPipelineReturnType } from 'store/slice/Pipeline/PipelineHook'
import { CurrentPipelineInfo } from 'components/common/CurrentPipelineInfo'

const FlowChart = React.memo<UseRunPipelineReturnType>((props) => {
  const open = useSelector(selectRightDrawerIsOpen)
  return (
    <RootDiv>
      <DndProvider backend={HTML5Backend}>
        <StyledDrawer variant="permanent">
          <MuiToolbar />
          <CurrentPipelineInfo />
          <DrawerContents>
            <AlgorithmTreeView />
          </DrawerContents>
        </StyledDrawer>

        <MainContents open={open}>
          <MuiToolbar />
          <ReactFlowComponent {...props} />
        </MainContents>
      </DndProvider>
      <RightDrawer />
    </RootDiv>
  )
})

export const drawerWidth = 240

const RootDiv = styled('div')({
  display: 'flex',
})

const StyledDrawer = styled(Drawer)({
  width: drawerWidth,
  flexShrink: 0,
  [`& .${drawerClasses.paper}`]: {
    width: drawerWidth,
  },
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
    marginRight: -rightDrawerWidth,
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
