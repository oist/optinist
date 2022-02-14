import React from 'react'
import { styled } from '@mui/material/styles'
import IconButton from '@mui/material/IconButton'
import Table from '@mui/material/Table'
import TableBody from '@mui/material/TableBody'
import TableCell from '@mui/material/TableCell'
import TableContainer from '@mui/material/TableContainer'
import TableHead from '@mui/material/TableHead'
import TableRow from '@mui/material/TableRow'
import Paper from '@mui/material/Paper'
import KeyboardArrowDownIcon from '@mui/icons-material/KeyboardArrowDown'
import KeyboardArrowUpIcon from '@mui/icons-material/KeyboardArrowUp'

import DoneIcon from '@mui/icons-material/Done'
import ErrorOutlineIcon from '@mui/icons-material/ErrorOutline'
import DeleteOutlineIcon from '@mui/icons-material/DeleteOutline'
import GetAppIcon from '@mui/icons-material/GetApp'

import { createData } from './DataType'
import { CollapsibleTable } from './CollapsibleTable'

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
      {/* <TableCell>Progress</TableCell> */}
      <TableCell>Import</TableCell>
      <TableCell>Delete</TableCell>
    </TableRow>
  )
}

const Row: React.FC<{
  row: ReturnType<typeof createData>
}> = ({ row }) => {
  const [open, setOpen] = React.useState(false)
  return (
    <React.Fragment>
      <TableRow
        sx={{
          '& > *': {
            borderBottom: 'unset',
          },
        }}
      >
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
        {/* <TableCell>
          <ProgressBar progress={row.progress} />
        </TableCell> */}
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

const ProgressBar: React.FC<{
  progress: number
}> = ({ progress }) => {
  const valueInPercent = progress
  return (
    <ProgressBarRoot>
      <ProgressValue>{`${valueInPercent.toLocaleString()} %`}</ProgressValue>
      <Bar valueInPercent={valueInPercent} />
    </ProgressBarRoot>
  )
}

const ProgressBarRoot = styled('div')(({ theme }) => ({
  border: `1px solid ${theme.palette.divider}`,
  position: 'relative',
  overflow: 'hidden',
  width: '100%',
  height: 26,
  borderRadius: 2,
}))

const ProgressValue = styled('div')({
  position: 'absolute',
  lineHeight: '24px',
  width: '100%',
  display: 'flex',
  justifyContent: 'center',
})

const Bar = styled('div')<{ valueInPercent: number }>(
  ({ valueInPercent, theme }) => ({
    height: '100%',
    maxWidth: `${valueInPercent}%`,
    backgroundColor:
      valueInPercent === 100
        ? theme.palette.success.light
        : valueInPercent >= 50 && valueInPercent < 100
        ? theme.palette.warning.light
        : theme.palette.error.light,
  }),
)
