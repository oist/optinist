import { FC } from "react"
import { useSelector } from "react-redux"

import Box from "@mui/material/Box"

import { VisualizeItem } from "components/Workspace/Visualize/VisualizeItem"
import { VisualizeItemAddButton } from "components/Workspace/Visualize/VisualizeItemAddButton"
import { selectVisualizeItemLayout } from "store/slice/VisualizeItem/VisualizeItemSelectors"
import { twoDimarrayEqualityFn } from "utils/EqualityUtils"

export const FlexItemList: FC = () => {
  const layout = useSelector(selectVisualizeItemLayout, twoDimarrayEqualityFn)
  return (
    <Box display="flex" flexWrap="wrap" flexDirection="column" p={1} m={1}>
      {layout.map((row, i) => (
        <Box key={i.toFixed()} display="flex" flexDirection="row">
          {row.map((itemId) => (
            <VisualizeItem itemId={itemId} key={itemId} />
          ))}
          <VisualizeItemAddButton itemId={row[row.length - 1]} />
        </Box>
      ))}
      <VisualizeItemAddButton />
    </Box>
  )
}
