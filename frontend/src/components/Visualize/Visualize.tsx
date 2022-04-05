import React from 'react'
import Drawer, { drawerClasses } from '@mui/material/Drawer'
import { default as MuiToolbar } from '@mui/material/Toolbar'
import { styled } from '@mui/material/styles'
import { drawerWidth } from 'components/FlowChart/FlowChart'
import { FlexItemList } from './VisualizeItems'
import { VisualizeItemEditor } from './VisualizeItemEditor'

const Visualize: React.FC = () => {
  return (
    <RootDiv>
      <StyledDrawer variant="permanent">
        <MuiToolbar />
        <StyledDrawerContents>
          <VisualizeItemEditor />
        </StyledDrawerContents>
      </StyledDrawer>
      <MainContents>
        <MuiToolbar />
        <FlexItemList />
      </MainContents>
    </RootDiv>
  )
}

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

const StyledDrawerContents = styled('div')({
  overflow: 'auto',
})

const MainContents = styled('main')({
  display: 'flex',
  flexDirection: 'column',
  flexGrow: 1,
  height: '100vh',
})

export default Visualize
