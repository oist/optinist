import React from 'react'
import { useSelector, useDispatch } from 'react-redux'
import { createStyles, Theme, makeStyles } from '@material-ui/core/styles'
import Drawer from '@material-ui/core/Drawer'
import Toolbar from '@material-ui/core/Toolbar'
import Typography from '@material-ui/core/Typography'
import Divider from '@material-ui/core/Divider'
import IconButton from '@material-ui/core/IconButton'
import ChevronRightIcon from '@material-ui/icons/ChevronRight'
import Box from '@material-ui/core/Box'
import {
  selectRightDrawerCurrentNodeId,
  selectRightDrawerIsOpen,
  selectRightDrawerMode,
} from 'store/slice/RightDrawer/RightDrawerSelectors'
import {
  closeRightDrawer,
  RIGHT_DRAWER_MODE,
} from 'store/slice/RightDrawer/RightDrawerSlice'
import { NWBSettingContents } from './NWB'
import { RootState } from 'store/store'
import { ParamForm } from './ParamForm'
import { ParamFormContext } from 'App'

const RightDrawer: React.FC = () => {
  const open = useSelector(selectRightDrawerIsOpen)
  const dispatch = useDispatch()
  const onClick = () => dispatch(closeRightDrawer())
  const classes = useStyles()
  const title = useSelector((state: RootState) => {
    const mode = selectRightDrawerMode(state)
    switch (mode) {
      case RIGHT_DRAWER_MODE.NWB:
        return 'NWB Setting'
      case RIGHT_DRAWER_MODE.PARAM_FORM:
        return 'Param From'
      default:
        return 'none'
    }
  })
  return (
    <Drawer
      open={open}
      anchor="right"
      variant="persistent"
      className={classes.drawer}
      classes={{
        paper: classes.drawerPaper,
      }}
    >
      <Toolbar />
      <Box display="flex" alignItems="center">
        <IconButton color="inherit" onClick={onClick}>
          <ChevronRightIcon />
        </IconButton>
        <Typography variant="h6">{title}</Typography>
      </Box>
      <Divider />
      <div className={classes.contents}>
        <Contents />
      </div>
    </Drawer>
  )
}

const Contents: React.FC = () => {
  const mode = useSelector(selectRightDrawerMode)
  switch (mode) {
    case RIGHT_DRAWER_MODE.NWB:
      return <NWBSettingContents />
    case RIGHT_DRAWER_MODE.PARAM_FORM:
      return <ParamFormConetent />
    default:
      return null
  }
}

const ParamFormConetent: React.FC = () => {
  const nodeId = useSelector(selectRightDrawerCurrentNodeId)
  if (nodeId != null) {
    return (
      <ParamFormContext.Provider value={nodeId}>
        <ParamForm />
      </ParamFormContext.Provider>
    )
  } else {
    return null
  }
}

export const rightDrawerWidth = 320

const useStyles = makeStyles(
  createStyles({
    drawer: {
      width: rightDrawerWidth,
      flexShrink: 0,
    },
    drawerPaper: {
      width: rightDrawerWidth,
    },
    contents: {
      height: '100%',
    },
  }),
)

export default RightDrawer
