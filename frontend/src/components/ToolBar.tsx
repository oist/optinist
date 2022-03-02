import React from 'react'
import Box from '@mui/material/Box'

import { UseRunPipelineReturnType } from 'store/slice/Pipeline/PipelineHook'
import { NWBSettingButton } from './FlowChart/NWB'
import { SnakemakeButton } from './FlowChart/Snakemake'
import { RunButtons } from './RunButtons'

export const ToolBar = React.memo<UseRunPipelineReturnType>((props) => (
  <Box
    display="flex"
    justifyContent="flex-end"
    sx={{
      width: '100%',
      p: 1,
      zIndex: (theme) => theme.zIndex.appBar,
    }}
  >
    <SnakemakeButton />
    <NWBSettingButton />
    <RunButtons {...props} />
  </Box>
))
