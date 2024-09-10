import { FC } from "react"
import { useSelector } from "react-redux"
import { useLocation } from "react-router-dom"

import MenuIcon from "@mui/icons-material/Menu"
import MuiAppBar from "@mui/material/AppBar"
import Box from "@mui/material/Box"
import IconButton from "@mui/material/IconButton"
import { styled } from "@mui/material/styles"
import Toolbar from "@mui/material/Toolbar"
import Typography from "@mui/material/Typography"

import Profile from "components/Layout/Profile"
import Tooltips from "components/Layout/Tooltips"
import Logo from "components/logo.png"
import WorkspaceTabs from "components/Workspace/WorkspaceTabs"
import { APP_BAR_HEIGHT } from "const/Layout"
import { selectModeStandalone } from "store/slice/Standalone/StandaloneSeclector"

const Header: FC<{
  handleDrawerOpen: () => void
}> = ({ handleDrawerOpen }) => {
  const isStandalone = useSelector(selectModeStandalone)
  return isStandalone ? (
    <StandaloneHeader />
  ) : (
    <MultiUserHeader handleDrawerOpen={handleDrawerOpen} />
  )
}

const StandaloneHeader: FC = () => {
  return (
    <StyledAppBar>
      <Toolbar>
        <Box display="flex" flexGrow={1}>
          <img src={Logo} alt="logo" width={75} height={50} />
        </Box>
        <WorkspaceTabs />
        <Tooltips />
      </Toolbar>
    </StyledAppBar>
  )
}

const MultiUserHeader: FC<{ handleDrawerOpen: () => void }> = ({
  handleDrawerOpen,
}) => {
  const showTabsRegex = /^\/console\/workspaces\/.+$/
  const location = useLocation()

  return (
    <StyledAppBar>
      <Toolbar>
        <IconButton
          color="inherit"
          aria-label="open drawer"
          onClick={handleDrawerOpen}
          edge="start"
        >
          <MenuIcon />
        </IconButton>
        <Box display="flex" flexGrow={1}>
          <TitleLogo>STUDIO</TitleLogo>
        </Box>
        {showTabsRegex.test(location.pathname) && <WorkspaceTabs />}
        <Profile />
        <Tooltips />
      </Toolbar>
    </StyledAppBar>
  )
}

const StyledAppBar = styled(MuiAppBar)({
  position: "fixed",
  backgroundColor: "#E1DEDB",
  color: "#000000",
  height: APP_BAR_HEIGHT,
})

const TitleLogo = styled(Typography)({
  fontWeight: 600,
  fontSize: 22,
})

export default Header
