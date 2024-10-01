import { DndProvider } from "react-dnd"
import { HTML5Backend } from "react-dnd-html5-backend"
import { TestBackend } from "react-dnd-test-backend"
import { Provider } from "react-redux"

import configureMockStore from "redux-mock-store"

import { describe, it } from "@jest/globals"
import { render } from "@testing-library/react"

import { ReactFlowComponent } from "components/Workspace/FlowChart/ReactFlowComponent"

const mockStore = configureMockStore()
const store = mockStore({
  flowElement: {
    flowNodes: [],
    flowEdges: [],
    elementCoord: { x: 0, y: 0 },
    loading: false,
  },
  pipeline: {
    run: {
      status: "StartUninitialized",
    },
  },
  currentPipeline: {
    uid: "test-uid",
  },
  runBtn: "RUN_NEW",
})

const isTestEnv = process.env.NODE_ENV === "test"

const backendFactory = isTestEnv ? TestBackend : HTML5Backend

describe("Flowchart Component", () => {
  // TODO: WIP - Fix the test
  // it("renders correctly", () => {
  //   const handleRunPipelineByUid = jest.fn()
  //   const handleRunPipeline = jest.fn()
  //   const handleCancelPipeline = jest.fn()
  //   render(
  //     <Provider store={store}>
  //       <DndProvider backend={backendFactory}>
  //         <ReactFlowComponent
  //           filePathIsUndefined={false}
  //           algorithmNodeNotExist={false}
  //           uid={undefined}
  //           status={"StartUninitialized"}
  //           runDisabled={false}
  //           handleRunPipelineByUid={handleRunPipelineByUid}
  //           handleRunPipeline={handleRunPipeline}
  //           handleCancelPipeline={handleCancelPipeline}
  //         />
  //       </DndProvider>
  //     </Provider>,
  //   )
  // })
})
