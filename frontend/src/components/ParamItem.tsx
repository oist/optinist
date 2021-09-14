import { useState, useContext, useEffect } from 'react'
import { makeStyles, Typography, Grid, Slider, Input } from '@material-ui/core'

import AppStateContext from 'contexts/AppStateContext'

const useStyles = makeStyles({
  root: {
    flexGrow: 1,
  },
  grid: {
    padding: 6,
  },
})

const ParamItem = (props: any) => {
  const default_value = 30
  const classes = useStyles()
  const [value, setValue] = useState<number | string>(default_value)
  const { state, dispatch } = useContext(AppStateContext)

  useEffect(() => {
    const currentAlgoParameters = state.algorithms.filter(function (algo) {
      return algo.name === state.currentSelectedAlgo
    })[0].parameters
    const currentParamValue = currentAlgoParameters.filter(function (param) {
      return param.name === props.name
    })[0].value

    setValue(currentParamValue)
  }, [state.currentSelectedAlgo])

  const handleSliderChange = (event: any, newValue: number | number[]) => {
    if (event.isTrusted) {
      setValue(newValue as number)
      dispatch({
        type: 'ParamUpdate',
        value: newValue as number,
        param: props.name,
      })
    } else {
      console.log(event)
    }
  }

  const handleInputChange = (event: React.ChangeEvent<HTMLInputElement>) => {
    setValue(event.target.value === '' ? '' : Number(event.target.value))
    dispatch({
      type: 'ParamUpdate',
      value: event.target.value === '' ? '' : Number(event.target.value),
      param: props.name,
    })
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
