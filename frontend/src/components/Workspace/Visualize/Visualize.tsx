import { FC, useState } from "react"

import { Box } from "@mui/material"
import { styled } from "@mui/material/styles"

import { CurrentPipelineInfo } from "components/common/CurrentPipelineInfo"
import { LeftSidebarContainer } from "components/common/LeftSidebarContainer"
import { FlexItemList } from "components/Workspace/Visualize/FlexItemList"
import { VisualizeProvider } from "components/Workspace/Visualize/VisualizeContext"
import { VisualizeItemEditor } from "components/Workspace/Visualize/VisualizeItemEditor"
import { CONTENT_HEIGHT } from "const/Layout"

const Visualize: FC = () => {
  const [leftSidebarOpen, setLeftSidebarOpen] = useState(true)

  return (
    <Box display="flex">
      <LeftSidebarContainer
        isOpen={leftSidebarOpen}
        onToggle={() => setLeftSidebarOpen(!leftSidebarOpen)}
      >
        <Box overflow="auto" marginRight={2}>
          <CurrentPipelineInfo />
          <VisualizeItemEditor />
        </Box>
      </LeftSidebarContainer>

      <MainContents>
        <VisualizeProvider>
          <FlexItemList />
        </VisualizeProvider>
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
