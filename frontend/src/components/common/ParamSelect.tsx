import { memo, ReactNode } from "react"

import {
  Box,
  FormControl,
  InputLabel,
  Select,
  SelectChangeEvent,
} from "@mui/material"

interface ParamSelectProps {
  label: string
  onChange: (event: SelectChangeEvent<string>) => void
  value?: string
  children: ReactNode
}

export const ParamSelect = memo(function ParamSelect({
  label,
  onChange,
  value,
  children,
}: ParamSelectProps) {
  return (
    <Box marginBottom={2}>
      <FormControl variant="standard" size="small" fullWidth>
        <InputLabel
          style={{
            position: "static",
            transform: "none",
            transition: "none",
            fontWeight: "bold",
            fontSize: "0.95rem",
            color: "text.secondary",
          }}
        >
          {label}
        </InputLabel>
        <Select variant="outlined" value={value} onChange={onChange}>
          {children}
        </Select>
      </FormControl>
    </Box>
  )
})
