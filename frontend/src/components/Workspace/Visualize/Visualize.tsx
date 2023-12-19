import { FC } from "react"

import { Box } from "@mui/material"
import { grey } from "@mui/material/colors"
import { styled } from "@mui/material/styles"

import { CurrentPipelineInfo } from "components/common/CurrentPipelineInfo"
import { FlexItemList } from "components/Workspace/Visualize/FlexItemList"
import { VisualizeItemEditor } from "components/Workspace/Visualize/VisualizeItemEditor"
import { CONTENT_HEIGHT, DRAWER_WIDTH } from "const/Layout"

const Visualize: FC = () => {
  return (
    <Box display="flex">
      <Box
        width={DRAWER_WIDTH}
        marginRight={3}
        borderRight={1}
        borderColor={grey[300]}
      >
        <Box overflow="auto" marginRight={2}>
          <CurrentPipelineInfo />
          <VisualizeItemEditor />
        </Box>
      </Box>

      <MainContents>
        <FlexItemList />
      </MainContents>
    </Box>
  )
}

const MainContents = styled("main")({
  display: "flex",
  flexDirection: "column",
  flexGrow: 1,
  minHeight: CONTENT_HEIGHT,
})

export default Visualize
