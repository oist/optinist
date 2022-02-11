import React from 'react'
import Box from '@material-ui/core/Box'
import Collapse from '@material-ui/core/Collapse'
import Table from '@material-ui/core/Table'
import TableBody from '@material-ui/core/TableBody'
import TableCell from '@material-ui/core/TableCell'
import TableHead from '@material-ui/core/TableHead'
import TableRow from '@material-ui/core/TableRow'
import Typography from '@material-ui/core/Typography'
import DoneIcon from '@material-ui/icons/Done'
import ErrorOutlineIcon from '@material-ui/icons/ErrorOutline'
import { createData } from './DataType'

export const CollapsibleTable: React.FC<{
  row: ReturnType<typeof createData>
  open: boolean
}> = ({ row, open }) => {
  return (
    <TableRow>
      <TableCell style={{ paddingBottom: 0, paddingTop: 0 }} colSpan={6}>
        <Collapse in={open} timeout="auto" unmountOnExit>
          <Box margin={1}>
            <Typography variant="h6" gutterBottom component="div">
              Details
            </Typography>
            <Table size="small" aria-label="purchases">
              <Head />
              <Body row={row} />
            </Table>
          </Box>
        </Collapse>
      </TableCell>
    </TableRow>
  )
}

const Head: React.FC = () => {
  return (
    <TableHead>
      <TableRow>
        <TableCell>Function</TableCell>
        <TableCell>Success</TableCell>
      </TableRow>
    </TableHead>
  )
}

const Body: React.FC<{
  row: ReturnType<typeof createData>
}> = ({ row }) => {
  return (
    <TableBody>
      {row.details.map((detailsRow) => (
        <TableRow key={detailsRow.function}>
          <TableCell component="th" scope="row">
            {detailsRow.function}
          </TableCell>
          <TableCell>
            {row.status ? (
              <DoneIcon style={{ color: 'green' }} />
            ) : (
              <ErrorOutlineIcon style={{ color: 'red' }} />
            )}
          </TableCell>
        </TableRow>
      ))}
    </TableBody>
  )
}
