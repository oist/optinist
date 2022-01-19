import React from 'react'
import PlotlyChart from 'react-plotlyjs-ts'
import { LinearProgress, Typography } from '@material-ui/core'
import { useSelector, useDispatch } from 'react-redux'

import { DisplayDataContext } from '../DisplayDataItem'
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
import { selectNodeLabelById } from 'store/slice/FlowElement/FlowElementSelectors'
import { DataGrid, GridRowsProp, GridColDef } from '@mui/x-data-grid'
import { RootState } from 'store/store'

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
  const { filePath: path, nodeId } = React.useContext(DisplayDataContext)
  // const label = useSelector(selectNodeLabelById(nodeId)
  const label = useSelector((state: RootState) => {
    if (nodeId) {
      return selectNodeLabelById(nodeId)(state)
    } else {
      return path
    }
  })
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
  // const sample_columns = useSelector(selectTableDataColumns(path), (a, b) => {
  //   if (a != null && b != null) {
  //     return arrayEqualityFn(a, b)
  //   } else {
  //     return a === undefined && b === undefined
  //   }
  // })

  // const rows: GridRowsProp = [
  //   { id: 1, col1: 'Hello', col2: 'World' },
  //   { id: 2, col1: 'DataGridPro', col2: 'is Awesome' },
  //   { id: 3, col1: 'Material-UI', col2: 'is Amazing' },
  // ];

  // const columns: GridColDef[] =
  //   useSelector(selectTableDataColumns(path), (a, b) => {
  //     if (a != null && b != null) {
  //       return arrayEqualityFn(a, b)
  //     } else {
  //       return a === undefined && b === undefined
  //     }
  //   }).map((x, idx) => {
  //     return {field: "col"+String(idx), headerName: x, width: 150}
  //   })

  // console.log(columns)
  // const sample_rows =
  //   tableData.map((row, row_id) => {
  //     // console.log(row)
  //     row.reduce((array, b, col_id) => {
  //       // console.log(array)
  //       return array
  //       // if (typeof array === "number") {
  //       //   return {id: row_id, array, ["col"+col_id]: b}
  //       // }
  //       // else if (typeof array !== "number") {
  //       //   return {array, ["col"+col_id]: b}
  //       // }
  //     });
  //   })

  //   tableData = [
  //     [1, 2, 3,],
  //     [4, 5, 6],
  //   ]

  //   [
  //     {id: 1, col1: 1, col2: 2, col3: 3,},
  //     {id: 2, col1: 1, col2: 2, col3: 3,},
  //   ]

  // console.log(rows)

  // const columns: GridColDef[] = [
  //   { field: 'col1', headerName: 'Column 1', width: 150 },
  //   { field: 'col2', headerName: 'Column 2', width: 150 },
  //   { field: 'col3', headerName: 'Column 1', width: 150 },
  // ];

  // return (
  //   <div style={{ height: 300, width: '100%' }}>
  //     <DataGrid rows={rows} columns={columns} />
  //   </div>
  // );

  const colmuns = useSelector(selectTableDataColumns(path), (a, b) => {
    if (a != null && b != null) {
      return arrayEqualityFn(a, b)
    } else {
      return a === undefined && b === undefined
    }
  })

  const data = React.useMemo(
    () =>
      colmuns === undefined || tableData === undefined
        ? []
        : [
            {
              type: 'table',
              header: {
                values: colmuns,
                align: 'center',
                fill: { color: 'gainsboro' },
              },
              cells: {
                values: transpose(tableData), // PlotlyChartの行列の向きに合わせるため転置する
              },
            },
          ],
    [colmuns, tableData],
  )
  const config = {
    displayModeBar: true,
  }
  const layout = {
    title: label,
    margin: {
      t: 30, // top
      l: 90, // left
      b: 30, // bottom
    },
    autosize: true,
    height: 350,
    yaxis: {
      autorange: 'reversed',
    },
  }
  return <PlotlyChart data={data} layout={layout} config={config} />
})

function transpose(array: number[][]) {
  return array[0].map((_, c) => array.map((r) => r[c]))
}
