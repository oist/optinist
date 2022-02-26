import React from 'react'
import { default as MuiToolbar } from '@mui/material/Toolbar'
import { ExperimentTable } from './ExperimentTable'

const Experiment: React.FC = () => {
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
        {/* <VisualizeItems /> */}
        <ExperimentTable />
      </main>
    </div>
  )
}

export default Experiment
