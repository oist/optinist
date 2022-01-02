import React from 'react'
import PlotlyChart from 'react-plotlyjs-ts'
import { LinearProgress, Typography } from '@material-ui/core'
import { useSelector, useDispatch } from 'react-redux'

import { DisplayDataTabContext } from 'App'
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

export const TablePlot = React.memo(() => {
  const { filePath: path } = React.useContext(DisplayDataTabContext)
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
  const { filePath: path, nodeId } = React.useContext(DisplayDataTabContext)
  const label = useSelector(selectNodeLabelById(nodeId))
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
