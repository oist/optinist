import React from 'react'
import IconButton from '@mui/material/IconButton'
import Close from '@mui/icons-material/Close'
import { SnackbarProvider, SnackbarKey, useSnackbar } from 'notistack'
import { BrowserRouter, Navigate, Route, Routes } from 'react-router-dom'
import Layout from 'components/Layout'
import Dashboard from 'pages/Dashboard'
import Account from 'pages/Account'
import AccountManager from 'pages/AccountManager'
import AccountDelete from 'pages/AccountDelete'
import Login from 'pages/Login'
import ResetPassword from 'pages/ResetPassword'
import Workspaces from 'pages/Workspace'
import Workspace from 'pages/Workspace/Workspace'
import { IS_STANDALONE } from 'const/Mode'

const App: React.FC = () => {
  return (
    <SnackbarProvider
      maxSnack={5}
      action={(snackbarKey) => (
        <SnackbarCloseButton snackbarKey={snackbarKey} />
      )}
    >
      <BrowserRouter>
        <Layout>
          {IS_STANDALONE ? (
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
              <Route path="/console/account-manager" element={<AccountManager />} />
              <Route path="/console/workspaces">
                <Route path="" element={<Workspaces />} />
                <Route path=":workspaceId" element={<Workspace />} />
              </Route>
              <Route path="/console/*" element={<Navigate replace to="/console" />} />
              <Route path="*" element={<Navigate replace to="/" />} />
            </Routes>
          )}
        </Layout>
      </BrowserRouter>
    </SnackbarProvider>
  )
}

const SnackbarCloseButton: React.FC<{ snackbarKey: SnackbarKey }> = ({
  snackbarKey,
}) => {
  const { closeSnackbar } = useSnackbar()
  return (
    <IconButton onClick={() => closeSnackbar(snackbarKey)} size="large">
      <Close style={{ color: 'white' }} />
    </IconButton>
  )
}

export default App
