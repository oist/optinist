import React from 'react'

import { useDispatch } from 'react-redux'

import Box from '@material-ui/core/Box'
import Paper from '@material-ui/core/Paper'
import AddIcon from '@material-ui/icons/Add'
import { makeStyles, createStyles, Theme } from '@material-ui/core/styles'
import Button from '@material-ui/core/Button'
import { addInitialItem } from 'store/slice/VisualizeItem/VisualizeItemSlice'

export const VisualizeItemAddButton: React.FC = () => {
  const dispatch = useDispatch()
  const onClick = () => {
    dispatch(addInitialItem())
  }
  const classes = useStyles()
  return (
    <Paper elevation={1} variant="outlined" className={classes.paper}>
      <Box
        display="flex"
        justifyContent="center"
        alignItems="center"
        height="100%"
      >
        <Button onClick={onClick} className={classes.button}>
          <AddIcon fontSize="large" color="primary" />
        </Button>
      </Box>
    </Paper>
  )
}

const useStyles = makeStyles((theme: Theme) =>
  createStyles({
    paper: {
      width: 260,
      height: 255,
      border: 'dashed',
      borderWidth: 2,
      borderColor: theme.palette.divider,
      margin: theme.spacing(1),
    },
    button: {
      width: '100%',
      height: '100%',
    },
  }),
)
