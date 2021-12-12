import { LinearProgress } from '@material-ui/core'
import { TableDataContext } from 'App'
import React from 'react'
import PlotlyChart from 'react-plotlyjs-ts'

import { useSelector, useDispatch } from 'react-redux'
import {
  csvIsUploadedByIdSelector,
  csvPathByIdSelector,
} from 'store/slice/FileData/FileDataSelector'
import { getTableData } from 'store/slice/PlotData/PlotDataAction'
import {
  tableDataColumnsSelector,
  tableDataIsLoadedSelector,
  tableDataSelector,
} from 'store/slice/PlotData/PlotDataSelector'
import { TableData } from 'store/slice/PlotData/PlotDataType'
import { RootState } from 'store/store'
import { arrayEqualityFn, twoDimarrayEqualityFn } from 'utils/EqualityUtils'

export const TablePlot = React.memo(() => {
  const { nodeId } = React.useContext(TableDataContext)
  const isUploaded = useSelector(csvIsUploadedByIdSelector(nodeId))
  if (isUploaded === true) {
    return <TablePlotContainer />
  } else {
    return null
  }
})

const TablePlotContainer = React.memo(() => {
  const { nodeId } = React.useContext(TableDataContext)
  const dispatch = useDispatch()
  const path = useSelector((state: RootState) => {
    return csvPathByIdSelector(nodeId)(state)
  })
  const isLoaded = useSelector(
    tableDataIsLoadedSelector(path ?? ''), // 応急処置
  )
  React.useEffect(() => {
    if (!isLoaded && path != null) {
      dispatch(getTableData({ path }))
    }
  }, [dispatch, isLoaded, path])
  if (isLoaded && path != null) {
    return <TablePlotImple path={path} />
  } else if (!isLoaded && path != null) {
    return <LinearProgress />
  } else {
    return null
  }
})

const TablePlotImple = React.memo<{ path: string }>(({ path }) => {
  const tableData = useSelector(
    tableDataSelector(path),
    (a: TableData | undefined, b: TableData | undefined) => {
      if (a != null && b != null) {
        return twoDimarrayEqualityFn(a, b)
      } else {
        return a === undefined && b === undefined
      }
    },
  )
  const colmuns = useSelector(tableDataColumnsSelector(path), (a, b) => {
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
