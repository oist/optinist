import { FC, useEffect, useState } from "react"
import { useDispatch, useSelector } from "react-redux"
import { BrowserRouter, Navigate, Route, Routes } from "react-router-dom"

import { isAxiosError } from "axios"
import { SnackbarProvider, SnackbarKey, useSnackbar } from "notistack"

import Close from "@mui/icons-material/Close"
import IconButton from "@mui/material/IconButton"

import BackendUnavailable from "components/common/BackendUnavailable"
import Loading from "components/common/Loading"
import Layout from "components/Layout"
import { RETRY_MAX_COUNT, RETRY_WAIT, RETRY_WAIT_LONG } from "const/Mode"
import Account from "pages/Account"
import AccountDelete from "pages/AccountDelete"
import AccountManager from "pages/AccountManager"
import Dashboard from "pages/Dashboard"
import Login from "pages/Login"
import ResetPassword from "pages/ResetPassword"
import Workspaces from "pages/Workspace"
import Workspace from "pages/Workspace/Workspace"
import { getModeStandalone } from "store/slice/Standalone/StandaloneActions"
import {
  selectLoading,
  selectModeStandalone,
} from "store/slice/Standalone/StandaloneSeclector"
import { AppDispatch } from "store/store"

/**
 * Returns whether the error from getModeStandalone should be retried.
 *
 * - Treat as "backend down" (retry): network/timeout errors and HTTP 5xx.
 * - Treat as terminal (no retry): HTTP 4xx. Non-axios errors are retried for safety.
 */
const isRetryableBackendError = (error: unknown): boolean => {
  if (!isAxiosError(error)) return true
  if (!error.response) return true
  const status = error.response.status
  return status >= 500 && status < 600
}

const App: FC = () => {
  const dispatch = useDispatch<AppDispatch>()
  const isStandalone = useSelector(selectModeStandalone)
  const loading = useSelector(selectLoading)
  // Show the "backend unavailable" screen instead of the normal layout.
  const [showBackendError, setShowBackendError] = useState(false)
  // false = background polling has stopped (4xx); user must reload manually.
  const [isRetrying, setIsRetrying] = useState(true)

  useEffect(() => {
    let cancelled = false
    let timerId: ReturnType<typeof setTimeout> | undefined
    let retryCount = 0

    const getMode = () => {
      dispatch(getModeStandalone())
        .unwrap()
        .then(() => {
          if (cancelled) return
          // Recovered: hide the error screen.
          setShowBackendError(false)
          setIsRetrying(true)
        })
        .catch((error: unknown) => {
          if (cancelled) return
          if (!isRetryableBackendError(error)) {
            // Terminal error (4xx): stop polling, require manual reload.
            setShowBackendError(true)
            setIsRetrying(false)
            return
          }
          retryCount += 1
          const reachedMax = retryCount >= RETRY_MAX_COUNT
          if (reachedMax) {
            setShowBackendError(true)
          }
          const wait = reachedMax ? RETRY_WAIT_LONG : RETRY_WAIT
          timerId = setTimeout(() => {
            getMode()
          }, wait)
        })
    }

    getMode()

    return () => {
      cancelled = true
      if (timerId !== undefined) clearTimeout(timerId)
    }
    //eslint-disable-next-line
  }, [])

  if (showBackendError) {
    return <BackendUnavailable isRetrying={isRetrying} />
  }

  return loading ? (
    <Loading loading={true} />
  ) : (
    <SnackbarProvider
      maxSnack={5}
      action={(snackbarKey) => (
        <SnackbarCloseButton snackbarKey={snackbarKey} />
      )}
    >
      <BrowserRouter>
        <Layout>
          {isStandalone ? (
            <Routes>
              <Route path="/" element={<Workspace />} />
              <Route path="*" element={<Navigate replace to="/" />} />
            </Routes>
          ) : (
            <Routes>
              <Route path="/" element={<Navigate replace to="/console" />} />
              <Route path="/account-deleted" element={<AccountDelete />} />
              <Route path="/login" element={<Login />} />
              <Route path="/reset-password" element={<ResetPassword />} />
              <Route path="/console" element={<Dashboard />} />
              <Route path="/console/account" element={<Account />} />
              <Route
                path="/console/account-manager"
                element={<AccountManager />}
              />
              <Route path="/console/workspaces">
                <Route path="" element={<Workspaces />} />
                <Route path=":workspaceId" element={<Workspace />} />
              </Route>
              <Route
                path="/console/*"
                element={<Navigate replace to="/console" />}
              />
              <Route path="*" element={<Navigate replace to="/" />} />
            </Routes>
          )}
        </Layout>
      </BrowserRouter>
    </SnackbarProvider>
  )
}

const SnackbarCloseButton: FC<{ snackbarKey: SnackbarKey }> = ({
  snackbarKey,
}) => {
  const { closeSnackbar } = useSnackbar()
  return (
    <IconButton onClick={() => closeSnackbar(snackbarKey)} size="large">
      <Close style={{ color: "white" }} />
    </IconButton>
  )
}

export default App
