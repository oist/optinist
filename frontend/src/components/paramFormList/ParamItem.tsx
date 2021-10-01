import { useDispatch, useSelector } from 'react-redux'
import { updateParam } from 'redux/slice/Element/Element'
import {
  currentElementSelector,
  paramValueSelector,
} from 'redux/slice/Element/ElementSelector'
import { makeStyles, Typography, Grid, Slider, Input } from '@material-ui/core'

const useStyles = makeStyles({
  root: {
    flexGrow: 1,
  },
  grid: {
    padding: 6,
  },
})

const ParamItem = (props: { name: string }) => {
  const classes = useStyles()
  const currentElement = useSelector(currentElementSelector)
  const dispatch = useDispatch()
  const value = useSelector(paramValueSelector(currentElement, props.name))

  const handleSliderChange = (event: any, newValue: number | number[]) => {
    if (event.isTrusted) {
      if (typeof newValue === 'number') {
        dispatch(updateParam({ name: props.name, newValue }))
      }
    }
  }

  const handleInputChange = (event: React.ChangeEvent<HTMLInputElement>) => {
    const newValue = event.target.value === '' ? '' : Number(event.target.value)
    if (typeof newValue === 'number') {
      dispatch(updateParam({ name: props.name, newValue }))
    }
  }

  const handleBlur = () => {
    if (value < 0) {
      dispatch(updateParam({ name: props.name, newValue: 0 }))
    } else if (value > 100) {
      dispatch(updateParam({ name: props.name, newValue: 100 }))
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
