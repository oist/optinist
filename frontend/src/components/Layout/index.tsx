import { FC, ReactNode, useEffect, useState } from "react"
import { useSelector, useDispatch } from "react-redux"
import { useLocation, useNavigate } from "react-router-dom"

import { Box } from "@mui/material"
import { styled } from "@mui/material/styles"

import Loading from "components/common/Loading"
import Header from "components/Layout/Header"
import LeftMenu from "components/Layout/LeftMenu"
import { APP_BAR_HEIGHT } from "const/Layout"
import { selectModeStandalone } from "store/slice/Standalone/StandaloneSeclector"
import { getMe } from "store/slice/User/UserActions"
import { selectCurrentUser } from "store/slice/User/UserSelector"
import { AppDispatch } from "store/store"
import { getToken } from "utils/auth/AuthUtils"

const authRequiredPathRegex = /^\/console\/?.*/

const Layout = ({ children }: { children?: ReactNode }) => {
  const user = useSelector(selectCurrentUser)
  const location = useLocation()
  const navigate = useNavigate()
  const dispatch = useDispatch<AppDispatch>()
  const isStandalone = useSelector(selectModeStandalone)

  console.log("isStandalone", isStandalone)

  const [loading, setLoadingAuth] = useState(
    !isStandalone && authRequiredPathRegex.test(location.pathname),
  )

  useEffect(() => {
    !isStandalone &&
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
    const isLogin = location.pathname === "/login"

    try {
      if (token) {
        await dispatch(getMe())
        if (isLogin) navigate("/console")
        return
      } else if (!isLogin) throw new Error("fail auth")
    } catch {
      navigate("/login", { replace: true })
    } finally {
      if (loading) setLoadingAuth(false)
    }
  }

  if (loading) return <Loading />

  return isStandalone || authRequiredPathRegex.test(location.pathname) ? (
    <AuthedLayout>{children}</AuthedLayout>
  ) : (
    <UnauthedLayout>{children}</UnauthedLayout>
  )
}

const AuthedLayout: FC<{ children: ReactNode }> = ({ children }) => {
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

const UnauthedLayout: FC<{ children: ReactNode }> = ({ children }) => {
  return (
    <LayoutWrapper>
      <ContentBodyWrapper>
        <ChildrenWrapper>{children}</ChildrenWrapper>
      </ContentBodyWrapper>
    </LayoutWrapper>
  )
}

const LayoutWrapper = styled(Box)({
  height: "100%",
  width: "100%",
})

const ContentBodyWrapper = styled(Box)(() => ({
  backgroundColor: "#ffffff",
  display: "flex",
  paddingTop: APP_BAR_HEIGHT,
  height: `calc(100% - ${APP_BAR_HEIGHT}px)`,
  paddingRight: 10,
  overflow: "auto",
}))

const ChildrenWrapper = styled("main", {
  shouldForwardProp: (prop) => prop !== "open",
})(({ theme }) => ({
  flexGrow: 1,
  padding: theme.spacing(3),
  transition: theme.transitions.create("margin", {
    easing: theme.transitions.easing.sharp,
    duration: theme.transitions.duration.leavingScreen,
  }),
}))

export default Layout
