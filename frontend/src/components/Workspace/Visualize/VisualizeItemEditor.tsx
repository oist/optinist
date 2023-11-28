import { createContext, FC, useContext } from "react"
import { useSelector } from "react-redux"

import { Box } from "@mui/material"

import { CsvItemEditor } from "components/Workspace/Visualize/Editor/CsvItemEditor"
import { HeatmapItemEditor } from "components/Workspace/Visualize/Editor/HeatmapItemEditor"
import { ImageItemEditor } from "components/Workspace/Visualize/Editor/ImageItemEditor"
import { RoiItemEditor } from "components/Workspace/Visualize/Editor/RoiItemEditor"
import { SaveFig } from "components/Workspace/Visualize/Editor/SaveFig"
import { TimeSeriesItemEditor } from "components/Workspace/Visualize/Editor/TimeSeriesItemEditor"
import {
  DATA_TYPE,
  DATA_TYPE_SET,
} from "store/slice/DisplayData/DisplayDataType"
import {
  selectSelectedVisualizeItemId,
  selectVisualizeDataType,
} from "store/slice/VisualizeItem/VisualizeItemSelectors"

export const VisualizeItemEditor = () => {
  const selectedItemId = useSelector(selectSelectedVisualizeItemId)
  return (
    <>
      {selectedItemId != null ? (
        <SelectedItemIdContext.Provider value={selectedItemId}>
          <DisplayDataItemEditor />
        </SelectedItemIdContext.Provider>
      ) : (
        "Please select item..."
      )}
    </>
  )
}

export const SelectedItemIdContext = createContext<number>(NaN)

const DisplayDataItemEditor: FC = () => {
  const itemId = useContext(SelectedItemIdContext)
  const dataType = useSelector(selectVisualizeDataType(itemId))
  return (
    <Box marginTop={2} marginRight={2}>
      <DisplayEditor dataType={dataType} />
    </Box>
  )
}

const DisplayEditor: FC<{
  dataType: DATA_TYPE | null
}> = ({ dataType }) => {
  /* 他のtypeのEditorも必要になったら追加する */
  switch (dataType) {
    case DATA_TYPE_SET.IMAGE:
      return <ImageItemEditor />
    case DATA_TYPE_SET.CSV:
      return <CsvItemEditor />
    case DATA_TYPE_SET.HEAT_MAP:
      return <HeatmapItemEditor />
    case DATA_TYPE_SET.TIME_SERIES:
      return <TimeSeriesItemEditor />
    case DATA_TYPE_SET.ROI:
      return <RoiItemEditor />
    case DATA_TYPE_SET.SCATTER:
    case DATA_TYPE_SET.BAR:
    case DATA_TYPE_SET.HISTOGRAM:
    case DATA_TYPE_SET.LINE:
    case DATA_TYPE_SET.PIE:
    case DATA_TYPE_SET.POLAR:
      return <SaveFig />
    case DATA_TYPE_SET.HTML:
      return <div>html editor</div>
    default:
      return null
  }
}
