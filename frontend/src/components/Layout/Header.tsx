import { FC } from 'react'
import { useLocation } from 'react-router-dom'
import { styled } from '@mui/material/styles'
import MuiAppBar from '@mui/material/AppBar'
import Box from '@mui/material/Box'
import Toolbar from '@mui/material/Toolbar'
import Typography from '@mui/material/Typography'
import IconButton from '@mui/material/IconButton'
import MenuIcon from '@mui/icons-material/Menu'
import Logo from 'components/logo.png'
import Tooltips from 'components/Layout/Tooltips'
import WorkspaceTabs from 'components/Workspace/WorkspaceTabs'
import { IS_STANDALONE } from 'const/Mode'
import Profile from './Profile'

const Header: FC<{
  handleDrawerOpen: () => void
}> = ({ handleDrawerOpen }) => {
  return IS_STANDALONE ? (
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
  const showTabsRegex = /^\/workspaces\/.+$/
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
      </Toolbar>
    </StyledAppBar>
  )
}

const StyledAppBar = styled(MuiAppBar)({
  position: 'fixed',
  backgroundColor: '#E1DEDB',
  color: '#000000',
})

const TitleLogo = styled(Typography)({
  fontWeight: 600,
  fontSize: 22,
})

export default Header
