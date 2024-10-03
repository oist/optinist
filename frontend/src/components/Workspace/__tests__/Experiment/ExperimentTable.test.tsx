/* eslint-disable no-undef */
import { Provider } from "react-redux"

import configureStore from "redux-mock-store"

import { describe, it, beforeEach } from "@jest/globals"
import { Store, AnyAction } from "@reduxjs/toolkit"
import {
  render,
  screen,
  fireEvent,
  waitFor,
  within,
} from "@testing-library/react"

import { renameExperiment } from "api/experiments/Experiments"
import { ReproduceButton } from "components/Workspace/Experiment/Button/ReproduceButton"
import {
  ExperimentTable,
  ExperimentUidContext,
} from "components/Workspace/Experiment/ExperimentTable"
import {
  deleteExperimentByUid,
  getExperiments,
} from "store/slice/Experiments/ExperimentsActions"
import { selectExperimentList } from "store/slice/Experiments/ExperimentsSelectors"
import { Experiments } from "store/slice/Experiments/ExperimentsType"
import { reset } from "store/slice/VisualizeItem/VisualizeItemSlice"
import { reproduceWorkflow } from "store/slice/Workflow/WorkflowActions"
import { RootState } from "store/store"

jest.mock("api/experiments/Experiments", () => ({
  renameExperiment: jest.fn(), // Mock the renameExperiment function
}))

// Mock the getExperiments action
jest.mock("store/slice/Experiments/ExperimentsActions", () => ({
  getExperiments: jest.fn(() => ({ type: "getExperiments" })), // Return plain object action
}))

// Mock the deleteExperimentByUid and clearCurrentPipeline actions
jest.mock("store/slice/Experiments/ExperimentsActions", () => ({
  deleteExperimentByUid: jest.fn(),
}))

jest.mock("store/slice/Experiments/ExperimentsActions", () => ({
  getExperiments: jest.fn(() => ({ type: "getExperiments" })), // Return plain object action
}))

jest.mock("store/slice/Pipeline/PipelineSlice", () => ({
  clearCurrentPipeline: jest.fn(),
}))

// Mock necessary actions and selectors
jest.mock("store/slice/Workflow/WorkflowActions", () => ({
  reproduceWorkflow: jest.fn(() => ({ unwrap: jest.fn() })),
}))

jest.mock("store/slice/VisualizeItem/VisualizeItemSlice", () => ({
  reset: jest.fn(),
}))

jest.mock("store/slice/Experiments/ExperimentsSelectors", () => ({
  ...jest.requireActual("store/slice/Experiments/ExperimentsSelectors"),
  selectExperimentsStatusIsUninitialized: jest.fn(),
  selectExperimentsStatusIsError: jest.fn(),
  selectExperimentList: jest.fn(),
}))

// Mock necessary actions and API
jest.mock("api/experiments/Experiments", () => ({
  renameExperiment: jest.fn(),
}))

const mockStore = configureStore<RootState, AnyAction>([])

describe("ExperimentTable", () => {
  let store: Store<RootState, AnyAction>

  beforeEach(() => {
    store = mockStore({
      experiments: {
        status: "fulfilled",
        experimentList: {
          1: {
            uid: "1",
            name: "Experiment 1",
            functions: {
              caiman_mc_p4f6nsi9bh: {
                name: "caiman_mc",
                nodeId: "caiman_mc_p4f6nsi9bh",
                status: "success",
                hasNWB: true,
              },
            },
            startedAt: "2023-09-17",
            hasNWB: true,
          },
          2: {
            uid: "2",
            name: "Experiment 2",
            functions: {
              input_0: {
                name: "sample_mouse2p_image 2",
                nodeId: "input_0",
                status: "error",
                hasNWB: false,
              },
            },
            startedAt: "2023-09-15",
            hasNWB: true,
          },
        },
      } as Experiments,
      algorithmList: {
        isLatest: false,
        tree: {},
      },
      algorithmNode: {},
      displayData: {
        timeSeries: {},
        heatMap: {},
        image: {},
        csv: {},
        roi: {},
        scatter: {},
        bar: {},
        html: {},
        histogram: {},
        line: {},
        pie: {},
        polar: {},
        loading: false,
        statusRoi: {
          temp_add_roi: [],
          temp_delete_roi: [],
          temp_merge_roi: [],
        },
        loadingStack: [],
        isEditRoiCommitting: false,
      },
      fileUploader: {},
      mode: {
        mode: false,
        loading: false,
      },
      flowElement: {
        flowNodes: [],
        flowEdges: [],
        flowPosition: [0, 0, 0],
        elementCoord: {
          x: 0,
          y: 0,
        },
      },
      inputNode: {},
      handleColor: {
        colorMap: { "#000000": "#000000" },
        nextKey: 0,
      },
      filesTree: {},
      nwb: {
        params: {},
      },
      rightDrawer: {
        open: false,
        mode: "nwb",
        currendNodeId: null,
      },
      visualaizeItem: {
        selectedItemId: null,
        items: {},
        layout: [],
      },
      snakemake: {
        params: {},
      },
      pipeline: {
        run: {
          uid: "",
          status: "Finished",
          runPostData: {
            name: "",
            nodeDict: {},
            edgeDict: {},
            nwbParam: {},
            snakemakeParam: {},
            forceRunList: [],
          },
          runResult: {},
        },
        runBtn: 1,
      },
      hdf5: {
        isLoading: false,
        tree: [],
      },
      matlab: {
        isLoading: false,
        tree: [],
      },
      workspace: {
        workspace: {
          items: [],
          total: 0,
          limit: 0,
          offset: 0,
        },
        currentWorkspace: {
          statusRoi: undefined,
          roiFilePath: undefined,
          workspaceId: 1,
          workspaceName: undefined,
          selectedTab: 0,
          ownerId: undefined,
        },
        loading: false,
      },
      user: {
        loading: false,
      },
    })
    ;(selectExperimentList as jest.Mock).mockReturnValue([
      {
        uid: "1",
        name: "Experiment 1",
        functions: {},
        startedAt: "2023-09-17",
        hasNWB: true,
      },
      {
        uid: "2",
        name: "Experiment 2",
        functions: {},
        startedAt: "2023-09-15",
        hasNWB: true,
      },
    ])
  })

  it("should render the experiment list", () => {
    render(
      <Provider store={store}>
        <ExperimentTable />
      </Provider>,
    )
    // Check that the experiment names are rendered
    expect(screen.getByText("Experiment 1")).toBeInTheDocument()
    expect(screen.getByText("Experiment 2")).toBeInTheDocument()
  })

  it("Checkbox per Row Should be checked", () => {
    render(
      <Provider store={store}>
        <ExperimentTable />
      </Provider>,
    )

    const checkbox = screen.getAllByRole("checkbox")

    fireEvent.click(checkbox[1])

    expect(checkbox[1]).toBeChecked()
  })

  it("Select All Checkbox is available", () => {
    render(
      <Provider store={store}>
        <ExperimentTable />
      </Provider>,
    )

    // Find the <span> with the data-testid
    const selectAllCheckboxSpan = screen.getByTestId("select-all-checkbox")

    // Get the actual checkbox input inside the span
    const selectAllCheckbox = within(selectAllCheckboxSpan).getByRole(
      "checkbox",
    )

    // Initially, all checkboxes should be unchecked
    const checkboxes = screen.getAllByRole("checkbox")
    checkboxes.forEach((checkbox) => {
      expect(checkbox).not.toBeChecked()
    })

    // Click the "Select All" checkbox
    fireEvent.click(selectAllCheckbox)

    // Now all checkboxes should be checked
    checkboxes.forEach((checkbox) => {
      expect(checkbox).toBeChecked()
    })
  })

  it("expands and collapses the row when clicking the arrow button", async () => {
    render(
      <Provider store={store}>
        <ExperimentTable />
      </Provider>,
    )

    // Locate the arrow button (expand/collapse)
    const expandButtonFirst = screen.getAllByRole("button", {
      name: /expand row/i,
    })[0] // Assuming we're dealing with the first row

    // Initially, the details should not be rendered
    expect(screen.queryByText("Details")).not.toBeInTheDocument()

    // Click to expand
    fireEvent.click(expandButtonFirst)

    // The details should now be visible
    await waitFor(() => {
      expect(screen.getByText("Details")).toBeInTheDocument()
    })

    const doneIconList = screen.getAllByTestId("DoneIcon")

    // Verify that the Function 1 details are displayed
    expect(screen.getByText("caiman_mc")).toBeInTheDocument()
    expect(screen.getByText("caiman_mc_p4f6nsi9bh")).toBeInTheDocument()
    expect(doneIconList[0]).toBeInTheDocument()

    // Click to collapse
    fireEvent.click(expandButtonFirst)

    // The details should no longer be visible
    await waitFor(() => {
      expect(screen.queryByText("Details")).not.toBeInTheDocument()
    })
  })

  it("enables the Delete button when a checkbox is checked", () => {
    render(
      <Provider store={store}>
        <ExperimentTable />
      </Provider>,
    )

    // Ensure the Delete button is initially disabled
    const deleteButton = screen.getByTestId("delete-selected-button")
    expect(deleteButton).toBeDisabled()

    // Click on the first checkbox to select an experiment
    const checkbox = screen.getAllByRole("checkbox")[1] // Assuming the first checkbox is for the first row
    fireEvent.click(checkbox)

    // Verify that the checkbox is checked
    expect(checkbox).toBeChecked()

    // Now the Delete button should be enabled
    expect(deleteButton).toBeEnabled()
  })

  it("should sort experiments by Timestamp in descending order by default", () => {
    render(
      <Provider store={store}>
        <ExperimentTable />
      </Provider>,
    )

    const rows = screen.getAllByRole("row")

    expect(rows[1]).toHaveTextContent("2023-09-17")
    expect(rows[3]).toHaveTextContent("2023-09-15")
  })

  it("should sort experiments by name when the name header is clicked", () => {
    render(
      <Provider store={store}>
        <ExperimentTable />
      </Provider>,
    )

    const nameHeader = screen.getByText("Name")
    fireEvent.click(nameHeader)

    const rows = screen.getAllByRole("row")
    expect(rows[1]).toHaveTextContent("Experiment 1")
    expect(rows[3]).toHaveTextContent("Experiment 2")
  })

  it("should toggle sorting order when the name header is clicked twice", () => {
    render(
      <Provider store={store}>
        <ExperimentTable />
      </Provider>,
    )

    const nameHeader = screen.getByText("Name")
    fireEvent.click(nameHeader)
    fireEvent.click(nameHeader)

    const rows = screen.getAllByRole("row")
    expect(rows[1]).toHaveTextContent("Experiment 2")
    expect(rows[3]).toHaveTextContent("Experiment 1")
  })

  it("should show the delete confirmation dialog and cancel the deletion", async () => {
    render(
      <Provider store={store}>
        <ExperimentTable />
      </Provider>,
    )

    // Find the delete button in the first row
    const deleteButton = screen.getAllByTestId("delete-button")[0]

    // // Click the delete button
    fireEvent.click(deleteButton)

    // Wait for the dialog to appear
    const dialog = await waitFor(() =>
      screen.getByRole("dialog", { name: /delete record\?/i }),
    )
    expect(dialog).toBeInTheDocument()

    // Find and click the "Cancel" button
    const cancelButton = screen.getByText(/cancel/i)
    fireEvent.click(cancelButton)

    // Check if the dialog has been closed
    await waitFor(() => {
      expect(dialog).not.toBeInTheDocument()
    })
  })

  it("should show the reproduce confirmation dialog and cancel the reproduce", async () => {
    render(
      <Provider store={store}>
        <ExperimentTable />
      </Provider>,
    )

    // Find the reproduce button in the first row
    const reproduceButton = screen.getAllByTestId("reproduce-button")[0]

    // Click the reproduce button
    fireEvent.click(reproduceButton)

    const dialog = await waitFor(
      () => screen.getByRole("dialog"), // Remove the name filter for troubleshooting
    )

    // Wait for the dialog to appear
    // const dialog = await waitFor(() =>
    //   screen.getByRole("dialog", { name: /reproduce workflow\?/i }),
    // )
    expect(dialog).toBeInTheDocument()

    // Find and click the "Cancel" button
    const cancelButton = screen.getByText(/cancel/i)
    fireEvent.click(cancelButton)

    // Check if the dialog has been closed
    await waitFor(() => {
      expect(dialog).not.toBeInTheDocument()
    })
  })

  // TODO: WIP - Fix the test
  it.skip("renders and handles reproduction", async () => {
    render(
      <Provider store={store}>
        <ExperimentUidContext.Provider value="1">
          <ReproduceButton />
        </ExperimentUidContext.Provider>
      </Provider>,
    )

    // Simulate clicking the reproduce button to open the dialog
    fireEvent.click(screen.getByTestId("reproduce-button"))

    // Assert the confirm dialog appears
    expect(screen.getByText("Reproduce workflow?")).toBeInTheDocument()

    // Simulate confirming the reproduction
    fireEvent.click(screen.getByText("reproduce"))

    // Wait for the dispatch actions
    await waitFor(() => {
      expect(reproduceWorkflow).toHaveBeenCalledWith({
        workspaceId: 1,
        uid: "1",
      })
    })

    // Assert reset action was dispatched
    expect(reset).toHaveBeenCalled()
  })

  // TODO: WIP - Fix the test
  it.skip("should delete a single experiment when confirmed in the dialog", async () => {
    render(
      <Provider store={store}>
        <ExperimentTable />
      </Provider>,
    )

    // Find the delete button for the first experiment in the table row
    const deleteButton = screen.getAllByTestId("delete-button")[0] // Assuming each row has a delete button with this data-testid

    // Click the delete button in the row
    fireEvent.click(deleteButton)

    // Wait for the confirmation dialog to appear
    const dialog = await waitFor(() =>
      screen.getByRole("dialog", { name: /delete record\?/i }),
    )
    expect(dialog).toBeInTheDocument()

    // Find and click the "Confirm" button in the dialog
    const confirmButton = screen.getByTestId("delete-confirm-button")
    fireEvent.click(confirmButton)

    // Ensure the deleteExperimentByUid action is called with the correct UID
    expect(deleteExperimentByUid).toHaveBeenCalledWith("1")
  })

  // TODO: WIP - Fix the test
  it.skip("should dispatch the getExperiments action when the Reload button is clicked", () => {
    render(
      <Provider store={store}>
        <ExperimentTable />
      </Provider>,
    )

    // Find the Reload button by its text or other attributes
    const reloadButton = screen.getByText(/reload/i)

    // Click the Reload button
    fireEvent.click(reloadButton)

    // Ensure the getExperiments action is dispatched
    expect(getExperiments).toHaveBeenCalledTimes(1)
  })

  it.skip("allows renaming an experiment", async () => {
    render(
      <Provider store={store}>
        <ExperimentTable />
      </Provider>,
    )

    // Assert that the original name is rendered
    expect(screen.getByText("Experiment 1")).toBeInTheDocument()

    // Simulate clicking the name to enable editing
    fireEvent.click(screen.getByText("Experiment 1"))

    // The input should now be present
    const input = screen.getByPlaceholderText("Name")
    expect(input).toBeInTheDocument()

    // Simulate changing the name
    fireEvent.change(input, { target: { value: "New Experiment Name" } })

    // Simulate input losing focus (triggering onBlurEdit)
    fireEvent.blur(input)

    // Wait for the renameExperiment API call
    await waitFor(() => {
      expect(renameExperiment).toHaveBeenCalledWith(
        1,
        "1",
        "New Experiment Name",
      )
    })

    // After renaming, the UI should update with the new name
    expect(screen.getByText("New Experiment Name")).toBeInTheDocument()
  })
})
