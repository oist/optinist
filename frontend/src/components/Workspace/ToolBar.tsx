import React from 'react'
import Box from '@mui/material/Box'

import { UseRunPipelineReturnType } from 'store/slice/Pipeline/PipelineHook'
import { NWBSettingButton } from './FlowChart/NWB'
import { SnakemakeButton } from './FlowChart/Snakemake'
import { RunButtons } from './RunButtons'
import { ImportWorkflowConfigButton } from './FlowChart/ImportWorkflowConfigButton'

export const ToolBar = React.memo<UseRunPipelineReturnType>((props) => (
  <Box
    style={{
      position: 'absolute',
      float: 'right',
      textAlign: 'right',
      top: -7,
      right: 10,
      zIndex: 4,
      textTransform: 'none',
    }}
  >
    <ImportWorkflowConfigButton />
    <SnakemakeButton />
    <NWBSettingButton />
    <RunButtons {...props} />
  </Box>
))
