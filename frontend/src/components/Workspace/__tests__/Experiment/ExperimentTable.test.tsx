import { Provider } from "react-redux"

import configureStore from "redux-mock-store"

import { describe, it, beforeEach } from "@jest/globals"
import { Store, AnyAction } from "@reduxjs/toolkit"
import { render, screen, fireEvent } from "@testing-library/react"

import { ExperimentTable } from "components/Workspace/Experiment/ExperimentTable"
import { selectExperimentList } from "store/slice/Experiments/ExperimentsSelectors"
import { Experiments } from "store/slice/Experiments/ExperimentsType"
import { RootState } from "store/store"

jest.mock("store/slice/Experiments/ExperimentsSelectors", () => ({
  selectExperimentsStatusIsUninitialized: jest.fn(),
  selectExperimentsStatusIsFulfilled: jest.fn(),
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
          workspaceId: undefined,
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

  it("should update the experiment list after renaming a record", () => {
    render(
      <Provider store={store}>
        <ExperimentTable />
      </Provider>,
    )

    const inputElement = screen.getByLabelText("Experiment 1")
    fireEvent.change(inputElement, { target: { value: "New Name" } })

    expect(selectExperimentList(store.getState())).toEqual([
      {
        uid: "1",
        name: "New Name",
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
})
