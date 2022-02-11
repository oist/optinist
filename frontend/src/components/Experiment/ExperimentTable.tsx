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
import { createData } from './DataType'
import { CollapsibleTable } from './CollapsibleTable'

const useRowStyles = makeStyles({
  root: {
    '& > *': {
      borderBottom: 'unset',
    },
  },
})

const rows = [
  createData('2022-02-02', true, '100%', 'name1'),
  createData('2022-02-03', false, '80%', 'name2'),
  createData('2022-02-04', true, '100%', 'name3'),
  createData('2022-02-05', false, '80%', 'name4'),
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
      <TableCell>Status</TableCell>
      <TableCell>Progress</TableCell>
      <TableCell>Name</TableCell>
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
        <TableCell>
          {row.status ? (
            <DoneIcon style={{ color: 'green' }} />
          ) : (
            <ErrorOutlineIcon style={{ color: 'red' }} />
          )}
        </TableCell>
        <TableCell>{row.progress}</TableCell>
        <TableCell>{row.name}</TableCell>
      </TableRow>
      <CollapsibleTable row={row} open={open} />
    </React.Fragment>
  )
}
