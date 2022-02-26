import React from 'react'
import DoneIcon from '@mui/icons-material/Done'
import ErrorOutlineIcon from '@mui/icons-material/ErrorOutline'
import HorizontalRuleIcon from '@mui/icons-material/HorizontalRule'

import { EXPERIMENTS_STATUS } from 'store/slice/Experiments/ExperimentsType'

export const ExperimentStatusIcon = React.memo<{ status: EXPERIMENTS_STATUS }>(
  ({ status }) => {
    switch (status) {
      case 'error':
        return <ErrorOutlineIcon color="error" />
      case 'success':
        return <DoneIcon color="success" />
      case 'running':
        return <HorizontalRuleIcon color="inherit" />
    }
  },
)
