import React from 'react'
import Box from '@mui/material/Box'

import { UseRunPipelineReturnType } from 'store/slice/Pipeline/PipelineHook'
import { NWBSettingButton } from './FlowChart/NWB'
import { SnakemakeButton } from './FlowChart/Snakemake'
import { RunButtons } from './RunButtons'

export const ToolBar = React.memo<UseRunPipelineReturnType>((props) => (
  <div style={{ width: '100%' }}>
    <Box display="flex" justifyContent="flex-end" style={{ padding: 4 }}>
      <SnakemakeButton />
      <NWBSettingButton />
      <RunButtons {...props} />
    </Box>
  </div>
))
