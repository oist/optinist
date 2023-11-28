import { memo } from "react"

import { Grid, Switch, Typography } from "@mui/material"

interface ParamSwitchProps {
  name: string
  value: boolean
  onChange: () => void
}

export const ParamSwitch = memo(function ParamSwitch({
  name,
  onChange,
  value,
}: ParamSwitchProps) {
  return (
    <Typography component="div">
      <Grid container component="label" alignItems="center">
        <Grid item xs={9}>
          <Typography
            fontSize="0.95rem"
            fontWeight="bold"
            color="text.secondary"
          >
            {name}
          </Typography>
        </Grid>
        <Grid item xs={3}>
          <Switch checked={value} onChange={onChange} />
        </Grid>
      </Grid>
    </Typography>
  )
})
