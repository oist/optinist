import React from 'react'
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
        <ExperimentTable />
      </main>
    </div>
  )
})

export default Experiment
