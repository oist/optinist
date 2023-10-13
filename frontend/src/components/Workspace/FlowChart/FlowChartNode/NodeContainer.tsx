import React from 'react'
import { alpha, useTheme } from '@mui/material/styles'

export const NodeContainer: React.FC<{
  children: React.ReactNode
  selected: boolean
}> = ({ children, selected }) => {
  const theme = useTheme()
  return (
    <div
      style={{
        height: '100%',
        width: '100%',
        background: selected
          ? alpha(theme.palette.primary.light, 0.1)
          : undefined,
      }}
    >
      {children}
    </div>
  )
}
