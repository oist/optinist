import React from 'react'
import { createStyles, makeStyles } from '@material-ui/core/styles'
import { default as MuiToolbar } from '@material-ui/core/Toolbar'
import { drawerWidth } from 'components/FlowChart/FlowChart'
import { ExperimentTable } from './ExperimentTable'
// import { VisualizeItems } from './VisualizeItems'
// import { VisualizeItemEditor } from './VisualizeItemEditor'

const Experiment: React.FC = () => {
  const classes = useStyles()
  return (
    <div className={classes.root}>
      {/* <Drawer
        className={classes.drawer}
        variant="permanent"
        classes={{
          paper: classes.drawerPaper,
        }}
      >
        <MuiToolbar />
        <div className={classes.drawerContainer}>
          <VisualizeItemEditor />
        </div>
      </Drawer> */}
      <main className={classes.content}>
        <MuiToolbar />
        {/* <VisualizeItems /> */}
        <ExperimentTable />
      </main>
    </div>
  )
}

const useStyles = makeStyles(
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
    },
  }),
)

export default Experiment
