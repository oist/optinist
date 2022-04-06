import React from 'react'
import { default as MuiToolbar } from '@mui/material/Toolbar'
import { ExperimentTable } from './ExperimentTable'

const Experiment = React.memo(() => {
  return (
    <div style={{ display: 'flex' }}>
      <main
        style={{
          display: 'flex',
          flexDirection: 'column',
          flexGrow: 1,
          height: '100vh',
          padding: 16,
        }}
      >
        <MuiToolbar />
        <ExperimentTable />
      </main>
    </div>
  )
})

export default Experiment
