import { useState } from 'react'
import { makeStyles, Typography, Grid, Slider, Input } from '@material-ui/core'

const useStyles = makeStyles({
  root: {
    flexGrow: 1,
  },
  grid: {
    padding: 6,
  },
})

const ParamItem = (props: any) => {
  console.log(props)
  const default_value = 30
  const classes = useStyles()
  const [value, setValue] = useState<number | string>(default_value)

  const handleSliderChange = (event: any, newValue: number | number[]) => {
    console.log(event)
    setValue(newValue as number)
  }

  const handleInputChange = (event: React.ChangeEvent<HTMLInputElement>) => {
    setValue(event.target.value === '' ? '' : Number(event.target.value))
  }

  const handleBlur = () => {
    if (value < 0) {
      setValue(0)
    } else if (value > 100) {
      setValue(100)
    }
  }

  return (
    <div className={classes.root}>
      <Grid container spacing={1} className={classes.grid}>
        <Grid item xs={3}>
          <Typography>{props.name}</Typography>
        </Grid>
        <Grid item xs={5}>
          <Slider
            value={typeof value === 'number' ? value : 0}
            onChange={handleSliderChange}
            aria-labelledby="continuous-slider"
          />
        </Grid>
        <Grid item xs={3}>
          <Input
            value={value}
            onChange={handleInputChange}
            margin="dense"
            onBlur={handleBlur}
            inputProps={{
              step: 10,
              min: 0,
              max: 100,
              type: 'number',
              'aria-labelledby': 'input-slider',
            }}
          />
        </Grid>
      </Grid>
    </div>
  )
}

export default ParamItem
