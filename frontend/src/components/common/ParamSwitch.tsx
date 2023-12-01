import { memo } from "react"

import { Grid, Switch, Typography } from "@mui/material"

interface ParamSwitchProps {
  label: string
  value: boolean
  onChange: () => void
}

export const ParamSwitch = memo(function ParamSwitch({
  label,
  onChange,
  value,
}: ParamSwitchProps) {
  return (
    <Grid
      container
      component="label"
      alignItems="center"
      justifyContent="space-between"
      marginBottom={2}
    >
      <Grid item>
        <Typography fontSize="0.95rem" fontWeight="bold" color="text.secondary">
          {label}
        </Typography>
      </Grid>
      <Grid item>
        <Switch checked={value} onChange={onChange} />
      </Grid>
    </Grid>
  )
})
