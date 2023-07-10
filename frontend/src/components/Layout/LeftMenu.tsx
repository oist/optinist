import { FC } from 'react'
import { useNavigate } from 'react-router-dom'
import Drawer from '@mui/material/Drawer'
import List from '@mui/material/List'
import ListItem from '@mui/material/ListItem'
import ListItemButton from '@mui/material/ListItemButton'
import ListItemIcon from '@mui/material/ListItemIcon'
import ListItemText from '@mui/material/ListItemText'
import HomeIcon from '@mui/icons-material/Home'
import SourceIcon from '@mui/icons-material/Source'
import { DRAWER_WIDTH } from 'const/Layout'
import { Box } from '@mui/material'

const LeftMenu: FC<{ open: boolean; handleDrawerClose: () => void }> = ({
  open,
  handleDrawerClose,
}) => {
  const navigate = useNavigate()

  const onClickDashboard = () => {
    handleDrawerClose()
    navigate('/')
  }

  const onClickWorkspaces = () => {
    handleDrawerClose()
    navigate('/workspaces')
  }

  return (
    <>
      <Drawer anchor="left" open={open} onClose={handleDrawerClose}>
        <Box sx={{ width: DRAWER_WIDTH }}>
          <List>
            <ListItem key="dashboard" disablePadding>
              <ListItemButton onClick={onClickDashboard}>
                <ListItemIcon>
                  <HomeIcon />
                </ListItemIcon>
                <ListItemText primary="Dashboard" />
              </ListItemButton>
            </ListItem>
            <ListItem key="workspaces" disablePadding>
              <ListItemButton onClick={onClickWorkspaces}>
                <ListItemIcon>
                  <SourceIcon />
                </ListItemIcon>
                <ListItemText primary="Workspaces" />
              </ListItemButton>
            </ListItem>
          </List>
        </Box>
      </Drawer>
    </>
  )
}

export default LeftMenu
