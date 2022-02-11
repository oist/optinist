import React from 'react'
import { makeStyles } from '@material-ui/core/styles'
import IconButton from '@material-ui/core/IconButton'
import Table from '@material-ui/core/Table'
import TableBody from '@material-ui/core/TableBody'
import TableCell from '@material-ui/core/TableCell'
import TableContainer from '@material-ui/core/TableContainer'
import TableHead from '@material-ui/core/TableHead'
import TableRow from '@material-ui/core/TableRow'
import Paper from '@material-ui/core/Paper'
import KeyboardArrowDownIcon from '@material-ui/icons/KeyboardArrowDown'
import KeyboardArrowUpIcon from '@material-ui/icons/KeyboardArrowUp'

import DoneIcon from '@material-ui/icons/Done'
import ErrorOutlineIcon from '@material-ui/icons/ErrorOutline'
import DeleteOutlineIcon from '@material-ui/icons/DeleteOutline'
import GetAppIcon from '@material-ui/icons/GetApp'

import { createData } from './DataType'
import { CollapsibleTable } from './CollapsibleTable'
import clsx from 'clsx'
import { createStyles, createTheme } from '@material-ui/core'

const useRowStyles = makeStyles({
  root: {
    '& > *': {
      borderBottom: 'unset',
    },
  },
})

const rows = [
  createData('2022-02-02', 'name1', true, 100),
  createData('2022-02-03', 'name2', false, 80),
  createData('2022-02-04', 'name3', true, 100),
  createData('2022-02-05', 'name4', false, 30),
]

export const ExperimentTable: React.FC = () => {
  return (
    <TableContainer component={Paper}>
      <Table aria-label="collapsible table">
        <TableHead>
          <Head />
        </TableHead>
        <TableBody>
          {rows.map((row, idx) => (
            <Row key={idx} row={row} />
          ))}
        </TableBody>
      </Table>
    </TableContainer>
  )
}

const Head: React.FC = () => {
  return (
    <TableRow>
      <TableCell />
      <TableCell>Timestamp</TableCell>
      <TableCell>Name</TableCell>
      <TableCell>Status</TableCell>
      <TableCell>Progress</TableCell>
      <TableCell>Import</TableCell>
      <TableCell>Delete</TableCell>
    </TableRow>
  )
}

const Row: React.FC<{
  row: ReturnType<typeof createData>
}> = ({ row }) => {
  const [open, setOpen] = React.useState(false)
  const classes = useRowStyles()

  return (
    <React.Fragment>
      <TableRow className={classes.root}>
        <TableCell>
          <IconButton
            aria-label="expand row"
            size="small"
            onClick={() => setOpen(!open)}
          >
            {open ? <KeyboardArrowUpIcon /> : <KeyboardArrowDownIcon />}
          </IconButton>
        </TableCell>
        <TableCell component="th" scope="row">
          {row.date}
        </TableCell>
        <TableCell>{row.name}</TableCell>
        <TableCell>
          {row.status ? (
            <DoneIcon style={{ color: 'green' }} />
          ) : (
            <ErrorOutlineIcon style={{ color: 'red' }} />
          )}
        </TableCell>
        <TableCell>
          <ProgressBar progress={row.progress} />
        </TableCell>
        <TableCell>
          <GetAppIcon style={{ color: 'blue' }} />
        </TableCell>
        <TableCell>
          <DeleteOutlineIcon style={{ color: 'red' }} />
        </TableCell>
      </TableRow>
      <CollapsibleTable row={row} open={open} />
    </React.Fragment>
  )
}

const defaultTheme = createTheme()
const useStyles = makeStyles(
  (theme) =>
    createStyles({
      root: {
        border: `1px solid ${theme.palette.divider}`,
        position: 'relative',
        overflow: 'hidden',
        width: '100%',
        height: 26,
        borderRadius: 2,
      },
      value: {
        position: 'absolute',
        lineHeight: '24px',
        width: '100%',
        display: 'flex',
        justifyContent: 'center',
      },
      bar: {
        height: '100%',
        '&.low': {
          backgroundColor: '#f44336',
        },
        '&.medium': {
          backgroundColor: '#efbb5aa3',
        },
        '&.high': {
          backgroundColor: '#088208a3',
        },
      },
    }),
  { defaultTheme },
)

const ProgressBar: React.FC<{
  progress: number
}> = ({ progress }) => {
  const valueInPercent = progress
  const classes = useStyles()

  return (
    <div className={classes.root}>
      <div
        className={classes.value}
      >{`${valueInPercent.toLocaleString()} %`}</div>
      <div
        className={clsx(classes.bar, {
          low: valueInPercent < 50,
          medium: valueInPercent >= 50 && valueInPercent < 100,
          high: valueInPercent === 100,
        })}
        style={{ maxWidth: `${valueInPercent}%` }}
      />
    </div>
  )
}
