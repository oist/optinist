import React from 'react'
import PlotlyChart from 'react-plotlyjs-ts'
import { LinearProgress, Typography } from '@material-ui/core'
import { useSelector, useDispatch } from 'react-redux'

import { DisplayDataContext } from '../DataContext'
import { arrayEqualityFn, twoDimarrayEqualityFn } from 'utils/EqualityUtils'
import {
  selectTableData,
  selectTableDataColumns,
  selectTableDataError,
  selectTableDataIsFulfilled,
  selectTableDataIsInitialized,
  selectTableDataIsPending,
} from 'store/slice/DisplayData/DisplayDataSelectors'
import { getTableData } from 'store/slice/DisplayData/DisplayDataActions'
import { TableData } from 'store/slice/DisplayData/DisplayDataType'
import { DataGrid, GridColDef } from '@mui/x-data-grid'

export const TablePlot = React.memo(() => {
  const { filePath: path } = React.useContext(DisplayDataContext)
  const isInitialized = useSelector(selectTableDataIsInitialized(path))
  const isPending = useSelector(selectTableDataIsPending(path))
  const isFulfilled = useSelector(selectTableDataIsFulfilled(path))
  const error = useSelector(selectTableDataError(path))
  const dispatch = useDispatch()
  React.useEffect(() => {
    if (!isInitialized) {
      dispatch(getTableData({ path }))
    }
  }, [dispatch, isInitialized, path])
  if (isPending) {
    return <LinearProgress />
  } else if (error != null) {
    return <Typography color="error">{error}</Typography>
  } else if (isFulfilled) {
    return <TablePlotImple />
  } else {
    return null
  }
})

const TablePlotImple = React.memo(() => {
  const { filePath: path } = React.useContext(DisplayDataContext)

  const tableData = useSelector(
    selectTableData(path),
    (a: TableData | undefined, b: TableData | undefined) => {
      if (a != null && b != null) {
        return twoDimarrayEqualityFn(a, b)
      } else {
        return a === undefined && b === undefined
      }
    },
  )

  const columns: GridColDef[] = useSelector(
    selectTableDataColumns(path),
    (a, b) => {
      if (a != null && b != null) {
        return arrayEqualityFn(a, b)
      } else {
        return a === undefined && b === undefined
      }
    },
  ).map((x, idx) => {
    return { field: 'col' + String(idx), headerName: x, width: 150 }
  })

  const rows = tableData.map((row, row_id) => {
    let rowObj = Object.fromEntries(
      row.map((value, index) => {
        return [`col${index}`, value]
      }),
    )
    rowObj['id'] = row_id
    return rowObj
  })

  return (
    <div style={{ height: 300, width: '100%' }}>
      <DataGrid rows={rows} columns={columns} />
    </div>
  )
})
