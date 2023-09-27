import { FC, ReactNode, useEffect, useState } from 'react'
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
import Loading from 'components/common/Loading'
import { APP_BAR_HEIGHT } from 'const/Layout'

const authRequiredPathRegex = /^\/console\/?.*/

const Layout = ({ children }: { children?: ReactNode }) => {
  const user = useSelector(selectCurrentUser)
  const location = useLocation()
  const navigate = useNavigate()
  const dispatch = useDispatch()

  const [loading, setLoadingAuth] = useState(
    !IS_STANDALONE && authRequiredPathRegex.test(location.pathname),
  )

  useEffect(() => {
    !IS_STANDALONE &&
      authRequiredPathRegex.test(location.pathname) &&
      checkAuth()
    // eslint-disable-next-line react-hooks/exhaustive-deps
  }, [location.pathname, user])

  const checkAuth = async () => {
    if (user) {
      if (loading) setLoadingAuth(false)
      return
    }
    const token = getToken()
    const isLogin = location.pathname === '/login'

    try {
      if (token) {
        await dispatch(getMe())
        if (isLogin) navigate('/console')
        return
      } else if (!isLogin) throw new Error('fail auth')
    } catch {
      navigate('/login', { replace: true })
    } finally {
      if (loading) setLoadingAuth(false)
    }
  }

  if (loading) return <Loading />

  return IS_STANDALONE || authRequiredPathRegex.test(location.pathname) ? (
    <AuthedLayout>{children}</AuthedLayout>
  ) : (
    <UnauthedLayout>{children}</UnauthedLayout>
  )
}

const AuthedLayout: FC = ({ children }) => {
  const [open, setOpen] = useState(false)
  const handleDrawerOpen = () => {
    setOpen(true)
  }

  const handleDrawerClose = () => {
    setOpen(false)
  }
  return (
    <LayoutWrapper>
      <Header handleDrawerOpen={handleDrawerOpen} />
      <ContentBodyWrapper>
        <LeftMenu open={open} handleDrawerClose={handleDrawerClose} />
        <ChildrenWrapper>{children}</ChildrenWrapper>
      </ContentBodyWrapper>
    </LayoutWrapper>
  )
}

const UnauthedLayout: FC = ({ children }) => {
  return (
    <LayoutWrapper>
      <ContentBodyWrapper>
        <ChildrenWrapper>{children}</ChildrenWrapper>
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
  paddingTop: APP_BAR_HEIGHT,
  height: `calc(100% - ${APP_BAR_HEIGHT}px)`,
  paddingRight: 10,
  overflow: 'auto',
}))

const ChildrenWrapper = styled('main', {
  shouldForwardProp: (prop) => prop !== 'open',
})<{}>(({ theme }) => ({
  flexGrow: 1,
  padding: theme.spacing(3),
  transition: theme.transitions.create('margin', {
    easing: theme.transitions.easing.sharp,
    duration: theme.transitions.duration.leavingScreen,
  }),
}))

export default Layout
