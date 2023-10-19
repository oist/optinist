import React from "react"

import { Box } from "@mui/material"
import { grey } from "@mui/material/colors"
import { styled } from "@mui/material/styles"

import { CurrentPipelineInfo } from "components/common/CurrentPipelineInfo"
import { FlexItemList } from "components/Workspace/Visualize/FlexItemList"
import { VisualizeItemEditor } from "components/Workspace/Visualize/VisualizeItemEditor"
import { CONTENT_HEIGHT, DRAWER_WIDTH } from "const/Layout"

const Visualize: React.FC = () => {
  return (
    <RootDiv>
      <Box
        sx={{
          width: DRAWER_WIDTH,
        }}
        borderRight={1}
        borderColor={grey[300]}
      >
        <CurrentPipelineInfo />
        <StyledDrawerContents>
          <VisualizeItemEditor />
        </StyledDrawerContents>
      </Box>

      <MainContents>
        <FlexItemList />
      </MainContents>
    </RootDiv>
  )
}

const RootDiv = styled("div")({
  display: "flex",
})

const StyledDrawerContents = styled("div")({
  overflow: "auto",
})

const MainContents = styled("main")({
  display: "flex",
  flexDirection: "column",
  flexGrow: 1,
  minHeight: CONTENT_HEIGHT,
})

export default Visualize
