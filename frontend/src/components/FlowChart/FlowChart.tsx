import React from 'react'
import clsx from 'clsx'
import { useSelector } from 'react-redux'
import Drawer from '@material-ui/core/Drawer'
import { default as MuiToolbar } from '@material-ui/core/Toolbar'
import { createStyles, Theme, makeStyles } from '@material-ui/core/styles'
import { AlgorithmTreeView } from './TreeView'
import { ReactFlowComponent } from './ReactFlowComponent'
import { ToolBar } from '../ToolBar'
import RightDrawer, { rightDrawerWidth } from './RightDrawer'
import { selectRightDrawerIsOpen } from 'store/slice/RightDrawer/RightDrawerSelectors'
import { UseRunPipelineReturnType } from 'store/slice/Pipeline/PipelineHook'

const FlowChart = React.memo<UseRunPipelineReturnType>((props) => {
  const classes = useStyles()
  const open = useSelector(selectRightDrawerIsOpen)
  return (
    <div className={classes.root}>
      <Drawer
        className={classes.drawer}
        variant="permanent"
        classes={{
          paper: classes.drawerPaper,
        }}
      >
        <MuiToolbar />
        <div className={classes.drawerContainer}>
          <AlgorithmTreeView />
        </div>
      </Drawer>
      <main
        className={clsx(classes.content, {
          [classes.contentShift]: open,
        })}
      >
        <MuiToolbar />
        <ToolBar {...props} />
        <ReactFlowComponent />
      </main>
      <RightDrawer />
    </div>
  )
})

export const drawerWidth = 240

const useStyles = makeStyles((theme: Theme) =>
  createStyles({
    root: {
      display: 'flex',
    },
    drawer: {
      width: drawerWidth,
      flexShrink: 0,
    },
    drawerPaper: {
      width: drawerWidth,
    },
    drawerContainer: {
      overflow: 'auto',
    },
    content: {
      display: 'flex',
      flexDirection: 'column',
      flexGrow: 1,
      height: '100vh',
      transition: theme.transitions.create('margin', {
        easing: theme.transitions.easing.sharp,
        duration: theme.transitions.duration.leavingScreen,
      }),
      marginRight: -rightDrawerWidth,
    },
    contentShift: {
      transition: theme.transitions.create('margin', {
        easing: theme.transitions.easing.easeOut,
        duration: theme.transitions.duration.enteringScreen,
      }),
      marginRight: 0,
    },
  }),
)

export default FlowChart
