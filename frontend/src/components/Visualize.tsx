import React from 'react'
import Drawer from '@material-ui/core/Drawer'
import { default as MuiToolbar } from '@material-ui/core/Toolbar'
import { useAppDrawerStyles } from './FlowChart/FlowChart'

const Visualize: React.FC = () => {
  const classes = useAppDrawerStyles()
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
        <div className={classes.drawerContainer}>設定とか</div>
      </Drawer>
      <main className={classes.content}>
        <MuiToolbar />
        ここにDisplayData
      </main>
    </div>
  )
}

export default Visualize
