import React from 'react'
import { useSelector, useDispatch } from 'react-redux'
import { styled } from '@mui/material/styles'
import Drawer, { drawerClasses } from '@mui/material/Drawer'
import Toolbar from '@mui/material/Toolbar'
import Typography from '@mui/material/Typography'
import Divider from '@mui/material/Divider'
import IconButton from '@mui/material/IconButton'
import ChevronRightIcon from '@mui/icons-material/ChevronRight'
import Box from '@mui/material/Box'
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
import { AlgorithmParamForm } from './AlgorithmParamForm'
import { SnakemakeContents } from './Snakemake'

const RightDrawer: React.FC = () => {
  const open = useSelector(selectRightDrawerIsOpen)
  const dispatch = useDispatch()
  const onClick = () => dispatch(closeRightDrawer())
  const title = useSelector((state: RootState) => {
    const mode = selectRightDrawerMode(state)
    switch (mode) {
      case RIGHT_DRAWER_MODE.NWB:
        return 'NWB Setting'
      case RIGHT_DRAWER_MODE.PARAM_FORM:
        return 'Param From'
      case RIGHT_DRAWER_MODE.SNAKEMAKE:
        return 'Snakemake'
      default:
        return 'none'
    }
  })
  return (
    <StyledDrawer open={open} anchor="right" variant="persistent">
      <Toolbar />
      <Box display="flex" alignItems="center">
        <IconButton color="inherit" onClick={onClick} size="large">
          <ChevronRightIcon />
        </IconButton>
        <Typography variant="h6">{title}</Typography>
      </Box>
      <Divider />
      <MainContents>
        <Contents />
      </MainContents>
    </StyledDrawer>
  )
}

const Contents: React.FC = () => {
  const mode = useSelector(selectRightDrawerMode)
  switch (mode) {
    case RIGHT_DRAWER_MODE.NWB:
      return <NWBSettingContents />
    case RIGHT_DRAWER_MODE.PARAM_FORM:
      return <ParamFormConetent />
    case RIGHT_DRAWER_MODE.SNAKEMAKE:
      return <SnakemakeContents />
    default:
      return null
  }
}

/**
 * nodeId
 */
export const ParamFormContext = React.createContext<string>('')

const ParamFormConetent: React.FC = () => {
  const nodeId = useSelector(selectRightDrawerCurrentNodeId)
  if (nodeId != null) {
    return (
      <ParamFormContext.Provider value={nodeId}>
        <AlgorithmParamForm />
      </ParamFormContext.Provider>
    )
  } else {
    return null
  }
}

export const rightDrawerWidth = 320

const StyledDrawer = styled(Drawer)({
  width: rightDrawerWidth,
  flexShrink: 0,
  [`& .${drawerClasses.paper}`]: {
    width: rightDrawerWidth,
  },
})

const MainContents = styled('main')({
  height: '100%',
})

export default RightDrawer
