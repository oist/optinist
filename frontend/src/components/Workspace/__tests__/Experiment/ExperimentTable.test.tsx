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

import { ExperimentTable } from "components/Workspace/Experiment/ExperimentTable"
import { selectExperimentList } from "store/slice/Experiments/ExperimentsSelectors"
import { Experiments } from "store/slice/Experiments/ExperimentsType"
import { RootState } from "store/store"

jest.mock("api/experiments/Experiments", () => ({
  renameExperiment: jest.fn(), // Mock the renameExperiment function
}))

jest.mock("store/slice/Experiments/ExperimentsSelectors", () => ({
  ...jest.requireActual("store/slice/Experiments/ExperimentsSelectors"),
  selectExperimentsStatusIsUninitialized: jest.fn(),
  selectExperimentsStatusIsError: jest.fn(),
  selectExperimentList: jest.fn(),
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
            startedAt: "2023-09-17",
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
          status: "StartSuccess",
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
        startedAt: "2023-09-17",
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
})
