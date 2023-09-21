import React from 'react'
import { ExperimentTable } from './ExperimentTable'
import { CONTENT_HEIGHT } from 'const/Layout'

const Experiment = React.memo(() => {
  return (
    <div style={{ display: 'flex' }}>
      <main
        style={{
          display: 'flex',
          flexDirection: 'column',
          flexGrow: 1,
          minHeight: CONTENT_HEIGHT,
        }}
      >
        <ExperimentTable />
      </main>
    </div>
  )
})

export default Experiment
