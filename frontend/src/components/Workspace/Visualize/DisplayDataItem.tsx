import { memo } from "react"
import { useSelector } from "react-redux"

import { DisplayDataContext } from "components/Workspace/Visualize/DataContext"
import { BarPlot } from "components/Workspace/Visualize/Plot/BarPlot"
import { CsvPlot } from "components/Workspace/Visualize/Plot/CsvPlot"
import { HeatMapPlot } from "components/Workspace/Visualize/Plot/HeatMapPlot"
import { HistogramPlot } from "components/Workspace/Visualize/Plot/HistogramPlot"
import { HTMLPlot } from "components/Workspace/Visualize/Plot/HTMLPlot"
import { ImagePlot } from "components/Workspace/Visualize/Plot/ImagePlot"
import { LinePlot } from "components/Workspace/Visualize/Plot/LinePlot"
import { MatlabPlot } from "components/Workspace/Visualize/Plot/MatlabPlot"
import { PiePlot } from "components/Workspace/Visualize/Plot/PiePlot"
import { PolarPlot } from "components/Workspace/Visualize/Plot/PolarPlot"
import { RoiPlot } from "components/Workspace/Visualize/Plot/RoiPlot"
import { ScatterPlot } from "components/Workspace/Visualize/Plot/ScatterPlot"
import { TimeSeriesPlot } from "components/Workspace/Visualize/Plot/TimeSeriesPlot"
import {
  DATA_TYPE,
  DATA_TYPE_SET,
} from "store/slice/DisplayData/DisplayDataType"
import {
  selectVisualizeDataFilePath,
  selectVisualizeDataNodeId,
  selectVisualizeDataType,
} from "store/slice/VisualizeItem/VisualizeItemSelectors"

interface DisplayDataContextType {
  itemId: number
}

export const DisplayDataItem = memo(function DisplayDataItem({
  itemId,
}: DisplayDataContextType) {
  const filePath = useSelector(selectVisualizeDataFilePath(itemId))
  const nodeId = useSelector(selectVisualizeDataNodeId(itemId))
  const dataType = useSelector(selectVisualizeDataType(itemId))
  if (filePath != null && dataType != null) {
    return (
      <DisplayDataContext.Provider
        value={{ nodeId, filePath, dataType, itemId }}
      >
        <DisplayPlot dataType={dataType} />
      </DisplayDataContext.Provider>
    )
  } else {
    return <div>Please select item correctly.</div>
  }
})

interface DataTypeProps {
  dataType: DATA_TYPE
}

const DisplayPlot = memo(function DisplayPlot({ dataType }: DataTypeProps) {
  switch (dataType) {
    case DATA_TYPE_SET.CSV:
      return <CsvPlot />
    case DATA_TYPE_SET.MATLAB:
      return <MatlabPlot />
    case DATA_TYPE_SET.TIME_SERIES:
      return <TimeSeriesPlot />
    case DATA_TYPE_SET.HEAT_MAP:
      return <HeatMapPlot />
    case DATA_TYPE_SET.IMAGE:
      return <ImagePlot />
    case DATA_TYPE_SET.ROI:
      return <RoiPlot />
    case DATA_TYPE_SET.SCATTER:
      return <ScatterPlot />
    case DATA_TYPE_SET.BAR:
      return <BarPlot />
    case DATA_TYPE_SET.HTML:
      return <HTMLPlot />
    case DATA_TYPE_SET.HISTOGRAM:
      return <HistogramPlot />
    case DATA_TYPE_SET.LINE:
      return <LinePlot />
    case DATA_TYPE_SET.PIE:
      return <PiePlot />
    case DATA_TYPE_SET.POLAR:
      return <PolarPlot />
    default:
      return null
  }
})
