import React from 'react'

export const RunPipeLineContext = React.createContext<{
  runPipeLine: any
  result: any
}>({ runPipeLine: null, result: null })
