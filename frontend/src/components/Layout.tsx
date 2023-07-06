import { FC, useEffect, useState } from 'react'
import { useSelector, useDispatch } from 'react-redux'
import { Link, useLocation, useNavigate } from 'react-router-dom'
import { Box, Typography } from '@mui/material'
import { styled } from '@mui/material/styles'
import { KeyboardBackspace } from '@mui/icons-material'
import HomeIcon from '@mui/icons-material/Home'
import SourceIcon from '@mui/icons-material/Source'
import Header from './Header'
import { getToken } from 'utils/auth/AuthUtils'
import { selectCurrentUser } from 'store/slice/User/UserSelector'
import { getMe } from 'store/slice/User/UserActions'

export const drawerWidth = 240

const ignorePaths = ['/login', '/account-delete', '/reset-password']
const loginPaths = ['/login', '/reset-password']

const Layout: FC = ({ children }) => {
  const user = useSelector(selectCurrentUser)
  const location = useLocation()
  const [width, setWidth] = useState(drawerWidth)
  const navigate = useNavigate()
  const dispatch = useDispatch()

  const onResize = () => {
    setWidth(width === drawerWidth ? 54 : drawerWidth)
  }

  useEffect(() => {
    checkAuth()
    // eslint-disable-next-line react-hooks/exhaustive-deps
  }, [location.pathname, user])

  const checkAuth = async () => {
    if (user) return
    const token = getToken()
    const isPageLogin = loginPaths.includes(window.location.pathname)

    try {
      if (token) {
        dispatch(getMe())
        if (isPageLogin) navigate('/')
        return
      } else if (!isPageLogin) throw new Error('fail auth')
    } catch {
      navigate('/login')
    }
  }

  return (
    <LayoutWrapper>
      {ignorePaths.includes(location.pathname) ? null : <Header />}
      <ContentBodyWrapper>
        {ignorePaths.includes(location.pathname) ? null : (
          <MenuLeft onResize={onResize} width={width} />
        )}
        <ChildrenWrapper
          style={{
            width: `calc(100% - ${
              ignorePaths.includes(location.pathname) ? 0 : width + 10
            }px)`,
            height: '100%',
            overflow: 'auto',
          }}
        >
          {children}
        </ChildrenWrapper>
      </ContentBodyWrapper>
    </LayoutWrapper>
  )
}

const MenuLeft: FC<{ onResize: () => void; width: number }> = ({
  onResize,
  width,
}) => {
  const { pathname } = useLocation()
  const isClose = width !== drawerWidth
  return (
    <MenuLeftWrapper style={{ width, minWidth: width }}>
      <BoxBack>
        <ButtonBack
          onClick={onResize}
          style={{ transform: `rotate(${width === drawerWidth ? 0 : 180}deg)` }}
        >
          <BoxDivider />
          <KeyboardBackspaceIcon />
        </ButtonBack>
      </BoxBack>
      <MenuList>
        <LinkWrapper to="/">
          <MenuItem isClose={isClose} active={pathname === '/'}>
            <HomeIcon />
            <TypographyMenu style={{ opacity: Number(width === drawerWidth) }}>
              Dashboard
            </TypographyMenu>
          </MenuItem>
        </LinkWrapper>
        <LinkWrapper to="/workspaces">
          <MenuItem isClose={isClose} active={pathname.includes('/workspaces')}>
            <SourceIcon />
            <TypographyMenu style={{ opacity: Number(width === drawerWidth) }}>
              Workspaces
            </TypographyMenu>
          </MenuItem>
        </LinkWrapper>
      </MenuList>
    </MenuLeftWrapper>
  )
}

const LinkWrapper = styled(Link)(() => ({
  textDecoration: 'none',
}))

const LayoutWrapper = styled(Box)({
  height: '100%',
  width: '100%',
})

const ContentBodyWrapper = styled(Box)(() => ({
  backgroundColor: '#ffffff',
  display: 'flex',
  paddingTop: 48,
  height: 'calc(100% - 48px)',
  paddingRight: 10,
  overflow: 'hidden',
}))

const ChildrenWrapper = styled(Box)(() => ({
  height: 'calc(100% - 10px)',
  display: 'flex',
  paddingTop: 10,
  paddingLeft: 10,
}))

const MenuLeftWrapper = styled(Box)({
  height: '100%',
  backgroundColor: '#283237',
  overflow: 'auto',
  transition: 'all 0.3s',
})

const BoxBack = styled(Box)({
  width: '100%',
  height: 54,
  display: 'flex',
  justifyContent: 'flex-end',
})

const ButtonBack = styled(Box)({
  height: 54,
  width: 54,
  display: 'flex',
  alignItems: 'center',
  justifyContent: 'center',
  cursor: 'pointer',
})

const BoxDivider = styled(Box)({
  height: 15,
  width: 1,
  backgroundColor: '#ffffff',
  marginRight: -2,
})

const KeyboardBackspaceIcon = styled(KeyboardBackspace)({
  color: '#ffffff',
  fontSize: 20,
})

const MenuList = styled('ul')({
  margin: 0,
  padding: 0,
})

const MenuItem = styled('li', {
  shouldForwardProp: (props) => props !== 'isClose' && props !== 'active',
})<{ isClose: boolean; active: boolean }>(({ isClose, active }) => ({
  padding: '0 15px',
  color: '#ffffff',
  listStyle: 'none',
  height: 38,
  minHeight: 38,
  display: 'flex',
  alignItems: 'center',
  gap: 10,
  width: 'calc(100% - 30px)',
  minWidth: 'max-content',
  transition: 'all 0.3s',
  cursor: 'pointer',
  backgroundColor: active ? 'rgba(255,255,255,0.4) !important' : 'transparent',
  '&:hover': {
    transform: isClose
      ? 'scale(1.05) translateX(2px)'
      : 'scale(1.05) translateX(10px)',
    backgroundColor: 'rgba(255,255,255,0.2)',
  },
}))

const TypographyMenu = styled(Typography)({
  lineHeight: '20px',
  marginTop: 4,
  fontWeight: 500,
  transition: 'opacity 0.3s',
})

export default Layout
