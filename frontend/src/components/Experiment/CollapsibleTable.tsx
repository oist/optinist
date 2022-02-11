import React from 'react'
import Box from '@mui/material/Box'
import Collapse from '@mui/material/Collapse'
import Table from '@mui/material/Table'
import TableBody from '@mui/material/TableBody'
import TableCell from '@mui/material/TableCell'
import TableHead from '@mui/material/TableHead'
import TableRow from '@mui/material/TableRow'
import Typography from '@mui/material/Typography'
import DoneIcon from '@mui/icons-material/Done'
import ErrorOutlineIcon from '@mui/icons-material/ErrorOutline'
import { createData } from './DataType'

export const CollapsibleTable: React.FC<{
  row: ReturnType<typeof createData>
  open: boolean
}> = ({ row, open }) => {
  return (
    <TableRow>
      <TableCell sx={{ paddingBottom: 0, paddingTop: 0 }} colSpan={6}>
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
