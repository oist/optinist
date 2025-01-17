/* eslint-disable no-undef */
import { Provider } from "react-redux"

import configureStore from "redux-mock-store"

import { describe, it } from "@jest/globals"
import { render, screen } from "@testing-library/react"

import Workspace from "pages/Workspace/Workspace"

// Mocking components imported in Workspace with displayName
const MockExperiment = () => <div>Experiment Component</div>
MockExperiment.displayName = "Experiment"

const MockFlowChart = () => <div>FlowChart Component</div>
MockFlowChart.displayName = "FlowChart"

const MockVisualize = () => <div>Visualize Component</div>
MockVisualize.displayName = "Visualize"

// Mocking components imported in Workspace with inline anonymous functions
jest.mock("components/Workspace/Experiment/Experiment", () => {
  const MockExperiment = () => <div>Experiment Component</div>
  MockExperiment.displayName = "Experiment"
  return MockExperiment
})

jest.mock("components/Workspace/FlowChart/FlowChart", () => {
  const MockFlowChart = () => <div>FlowChart Component</div>
  MockFlowChart.displayName = "FlowChart"
  return MockFlowChart
})

jest.mock("components/Workspace/Visualize/Visualize", () => {
  const MockVisualize = () => <div>Visualize Component</div>
  MockVisualize.displayName = "Visualize"
  return MockVisualize
})

jest.mock("store/slice/Pipeline/PipelineHook", () => ({
  useRunPipeline: () => ({}),
}))

const mockStore = configureStore([])

describe("Workspace Component", () => {
  it("should render the correct content when each tab is active", () => {
    const initialState = {
      workspace: {
        currentWorkspace: {
          selectedTab: 0,
        },
      },
    }
    let store = mockStore(initialState)

    const { rerender } = render(
      <Provider store={store}>
        <Workspace />
      </Provider>,
    )

    // Verify that the FlowChart component is displayed when activeTab is 0
    expect(screen.getByText("FlowChart Component")).toBeInTheDocument()
    expect(screen.queryByText("Visualize Component")).not.toBeInTheDocument()
    expect(screen.queryByText("Experiment Component")).not.toBeInTheDocument()

    // Change activeTab to 1 and rerender
    store = mockStore({
      workspace: {
        currentWorkspace: {
          selectedTab: 1,
        },
      },
    })
    rerender(
      <Provider store={store}>
        <Workspace />
      </Provider>,
    )

    // Verify that the Visualize component is displayed when activeTab is 1
    expect(screen.getByText("Visualize Component")).toBeInTheDocument()
    expect(screen.queryByText("FlowChart Component")).not.toBeInTheDocument()
    expect(screen.queryByText("Experiment Component")).not.toBeInTheDocument()

    // Change activeTab to 2 and rerender
    store = mockStore({
      workspace: {
        currentWorkspace: {
          selectedTab: 2,
        },
      },
    })
    rerender(
      <Provider store={store}>
        <Workspace />
      </Provider>,
    )

    // Verify that the Experiment component is displayed when activeTab is 2
    expect(screen.getByText("Experiment Component")).toBeInTheDocument()
    expect(screen.queryByText("FlowChart Component")).not.toBeInTheDocument()
    expect(screen.queryByText("Visualize Component")).not.toBeInTheDocument()
  })
})
