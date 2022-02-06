import React from 'react'

import { Box, IconButton } from '@material-ui/core'
import Close from '@material-ui/icons/Close'
import { SnackbarProvider, SnackbarKey, useSnackbar } from 'notistack'

import { NWBSettingButton } from './FlowChart/NWB'

import { SnakemakeButton } from './FlowChart/Snakemake'
import { RunButtons } from './RunButtons'
import { UseRunPipelineReturnType } from 'store/slice/Pipeline/PipelineHook'

export const ToolBar = React.memo<UseRunPipelineReturnType>((props) => (
  <SnackbarProvider
    maxSnack={5}
    action={(snackbarKey) => <SnackbarCloseButton snackbarKey={snackbarKey} />}
  >
    <ToolBarImple {...props} />
  </SnackbarProvider>
))

const SnackbarCloseButton: React.FC<{ snackbarKey: SnackbarKey }> = ({
  snackbarKey,
}) => {
  const { closeSnackbar } = useSnackbar()
  return (
    <IconButton onClick={() => closeSnackbar(snackbarKey)}>
      <Close style={{ color: 'white' }} />
    </IconButton>
  )
}

const ToolBarImple = React.memo<UseRunPipelineReturnType>((props) => {
  return (
    <div style={{ width: '100%' }}>
      <Box display="flex" justifyContent="flex-end" style={{ padding: 4 }}>
        <SnakemakeButton />
        <NWBSettingButton />
        <RunButtons {...props} />
      </Box>
    </div>
  )
})
