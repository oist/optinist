import { Provider } from "react-redux"

import configureStore from "redux-mock-store"

import { describe, it, beforeEach } from "@jest/globals"
import { Store, AnyAction } from "@reduxjs/toolkit"
import {
  render,
  screen,
  fireEvent,
  prettyDOM,
  waitFor,
  within,
} from "@testing-library/react"

import { renameExperiment } from "api/experiments/Experiments"
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
            functions: {},
            startedAt: "2023-09-17",
            hasNWB: true,
          },
          2: {
            uid: "2",
            name: "Experiment 2",
            functions: {},
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
})
