import React from 'react'
import Box from '@mui/material/Box'

import { UseRunPipelineReturnType } from 'store/slice/Pipeline/PipelineHook'
import { NWBSettingButton } from './FlowChart/NWB'
import { SnakemakeButton } from './FlowChart/Snakemake'
import { RunButtons } from './RunButtons'
import { Button } from '@mui/material'
import ArrowBackIosIcon from '@mui/icons-material/ArrowBackIos'
import { useNavigate } from 'react-router-dom'
import { IS_STANDALONE } from 'const/Mode'
import { ImportWorkflowConfigButton } from './FlowChart/ImportWorkflowConfigButton'

export const ToolBar = React.memo<UseRunPipelineReturnType>((props) => {
  const navigate = useNavigate()
  return (
    <Box
      style={{
        position: 'absolute',
        float: 'right',
        textAlign: 'right',
        top: -7,
        right: 10,
        zIndex: 4,
        textTransform: 'none',
        fontSize: '1rem',
      }}
    >
      { !IS_STANDALONE &&
        (
          <Button onClick={() => navigate('/console/workspaces')}>
            <ArrowBackIosIcon />
            Workspaces
          </Button>
        )
      }
      <ImportWorkflowConfigButton />
      <SnakemakeButton />
      <NWBSettingButton />
      <RunButtons {...props} />
    </Box>
  )
})
