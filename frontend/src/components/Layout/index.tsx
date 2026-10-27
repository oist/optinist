import { FC, ReactNode, useEffect, useState } from "react"
import { useSelector, useDispatch } from "react-redux"
import { useLocation, useNavigate } from "react-router-dom"

import { Box } from "@mui/material"
import { styled } from "@mui/material/styles"

import Loading from "components/common/Loading"
import { LogsFloatingButton } from "components/common/LogsFloatingButton"
import Header from "components/Layout/Header"
import LeftMenu from "components/Layout/LeftMenu"
import ModalLogs from "components/Workspace/FlowChart/ModalLogs"
import { APP_BAR_HEIGHT } from "const/Layout"
import { selectLogsModalIsOpen } from "store/slice/LogsModal/LogsModalSelectors"
import { closeLogsModal } from "store/slice/LogsModal/LogsModalSlice"
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

  const [loading, setLoading] = useState(
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
      if (loading) setLoading(false)
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
      if (loading) setLoading(false)
    }
  }

  return isStandalone || authRequiredPathRegex.test(location.pathname) ? (
    <AuthedLayout>{children}</AuthedLayout>
  ) : (
    <>
      <Loading loading={loading} />
      <UnauthedLayout>{children}</UnauthedLayout>
    </>
  )
}

const AuthedLayout: FC<{ children: ReactNode }> = ({ children }) => {
  const dispatch = useDispatch<AppDispatch>()
  const [open, setOpen] = useState(false)
  const logsModalOpen = useSelector(selectLogsModalIsOpen)
  const isStandalone = useSelector(selectModeStandalone)

  const handleDrawerOpen = () => {
    setOpen(true)
  }

  const handleDrawerClose = () => {
    setOpen(false)
  }

  const handleLogsModalClose = () => {
    dispatch(closeLogsModal())
  }

  return (
    <LayoutWrapper>
      <Header handleDrawerOpen={handleDrawerOpen} />
      <ContentBodyWrapper>
        <LeftMenu open={open} handleDrawerClose={handleDrawerClose} />
        <ChildrenWrapper>{children}</ChildrenWrapper>
      </ContentBodyWrapper>
      {!isStandalone && <LogsFloatingButton />}
      {logsModalOpen && <ModalLogs isOpen onClose={handleLogsModalClose} />}
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
