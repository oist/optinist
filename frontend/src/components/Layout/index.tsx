import { FC, useEffect, useState } from 'react'
import { useSelector, useDispatch } from 'react-redux'
import { useLocation, useNavigate } from 'react-router-dom'
import { Box } from '@mui/material'
import { styled } from '@mui/material/styles'
import { getToken } from 'utils/auth/AuthUtils'
import { selectCurrentUser } from 'store/slice/User/UserSelector'
import { getMe } from 'store/slice/User/UserActions'
import Header from './Header'
import LeftMenu from './LeftMenu'
import { IS_STANDALONE } from 'const/Mode'

const ignorePaths = ['/login', '/account-delete', '/reset-password']
const loginPaths = ['/login', '/reset-password']

const Layout: FC = ({ children }) => {
  const user = useSelector(selectCurrentUser)
  const location = useLocation()
  const [open, setOpen] = useState(false)
  const navigate = useNavigate()
  const dispatch = useDispatch()

  const handleDrawerOpen = () => {
    setOpen(true)
  }

  const handleDrawerClose = () => {
    setOpen(false)
  }

  useEffect(() => {
    !IS_STANDALONE && checkAuth()
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
      {ignorePaths.includes(location.pathname) ? null : (
        <Header handleDrawerOpen={handleDrawerOpen} />
      )}
      <ContentBodyWrapper>
        {ignorePaths.includes(location.pathname) ? null : (
          <LeftMenu open={open} handleDrawerClose={handleDrawerClose} />
        )}
        <ChildrenWrapper open={open}>
          {children}
        </ChildrenWrapper>
      </ContentBodyWrapper>
    </LayoutWrapper>
  )
}

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

const ChildrenWrapper = styled('main', { shouldForwardProp: (prop) => prop !== 'open' })<{
  open?: boolean;
}>(({ theme, open }) => ({
  flexGrow: 1,
  padding: theme.spacing(3),
  transition: theme.transitions.create('margin', {
    easing: theme.transitions.easing.sharp,
    duration: theme.transitions.duration.leavingScreen,
  }),
  ...(open && {
    transition: theme.transitions.create('margin', {
      easing: theme.transitions.easing.easeOut,
      duration: theme.transitions.duration.enteringScreen,
    }),
    marginLeft: 0,
  }),
}));

export default Layout
