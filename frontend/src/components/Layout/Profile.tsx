import { FC, useState, MouseEvent } from "react"
import { useDispatch } from "react-redux"
import { useNavigate } from "react-router-dom"

import AccountCircleIcon from "@mui/icons-material/AccountCircle"
import Logout from "@mui/icons-material/Logout"
import PortraitIcon from "@mui/icons-material/Portrait"
import {
  ListItemIcon,
  ListItemText,
  Menu,
  MenuItem,
  Tooltip,
} from "@mui/material"
import IconButton from "@mui/material/IconButton"

import { logout } from "store/slice/User/UserSlice"

const Profile: FC = () => {
  const [anchorEl, setAnchorEl] = useState<null | HTMLElement>(null)
  const navigate = useNavigate()
  const dispatch = useDispatch()
  const handleMenu = (event: MouseEvent<HTMLElement>) => {
    setAnchorEl(event.currentTarget)
  }

  const handleCloseMenu = () => {
    setAnchorEl(null)
  }

  const onClickLogout = () => {
    setAnchorEl(null)
    dispatch(logout())
    navigate("/login")
  }

  const onClickAccount = () => {
    setAnchorEl(null)
    navigate("/console/account")
  }

  return (
    <>
      <Tooltip title="Profile">
        <IconButton
          aria-label="open profile menu"
          aria-haspopup="true"
          onClick={handleMenu}
        >
          <AccountCircleIcon />
        </IconButton>
      </Tooltip>
      <Menu
        id="profile-menu"
        anchorEl={anchorEl}
        anchorOrigin={{
          vertical: "bottom",
          horizontal: "right",
        }}
        keepMounted
        transformOrigin={{
          vertical: "top",
          horizontal: "right",
        }}
        open={Boolean(anchorEl)}
        onClose={handleCloseMenu}
      >
        <MenuItem onClick={onClickAccount}>
          <ListItemIcon>
            <PortraitIcon />
          </ListItemIcon>
          <ListItemText>Account Profile</ListItemText>
        </MenuItem>
        <MenuItem onClick={onClickLogout}>
          <ListItemIcon>
            <Logout />
          </ListItemIcon>
          <ListItemText>Sign Out</ListItemText>
        </MenuItem>
      </Menu>
    </>
  )
}

export default Profile
