import React from "react"
import { Provider } from "react-redux"

import configureStore from "redux-mock-store"
import thunk from "redux-thunk"

// import { expect, describe, it } from "@jest/globals"
import { render, screen, fireEvent } from "@testing-library/react"

import { AlgorithmTreeView } from "components/Workspace/FlowChart/TreeView"
import { getAlgoList } from "store/slice/AlgorithmList/AlgorithmListActions"

jest.mock("store/slice/AlgorithmList/AlgorithmListActions", () => ({
  getAlgoList: jest.fn(),
}))

const mockStore = configureStore([thunk])

describe("AlgorithmTreeView", () => {
  let store: ReturnType<typeof mockStore>

  beforeEach(() => {
    store = mockStore({
      AlgorithmList: {
        isLated: false,
        tree: {
          Algorithm: {
            children: {
              node1: {
                functionPath: "path/to/node1",
              },
              node2: {
                functionPath: "path/to/node2",
              },
            },
          },
        },
      },
      Pipeline: {
        latestUid: "workflow123",
      },
    })

    store.dispatch = jest.fn()
  })

  it("renders the AlgorithmTreeView component", () => {
    render(
      <Provider store={store}>
        <AlgorithmTreeView />
      </Provider>,
    )

    expect(screen.getByText("Data")).toBeInTheDocument()
    expect(screen.getByText("Algorithm")).toBeInTheDocument()
    expect(screen.getByText("node1")).toBeInTheDocument()
    expect(screen.getByText("node2")).toBeInTheDocument()
  })

  it("dispatches getAlgoList action when component mounts", () => {
    render(
      <Provider store={store}>
        <AlgorithmTreeView />
      </Provider>,
    )

    expect(store.dispatch).toHaveBeenCalledWith(getAlgoList())
  })

  it("adds a new algorithm node when clicking the add button", () => {
    render(
      <Provider store={store}>
        <AlgorithmTreeView />
      </Provider>,
    )

    const addButton = screen.getAllByLabelText("add")[0]
    fireEvent.click(addButton)

    expect(store.dispatch).toHaveBeenCalledWith(
      expect.objectContaining({
        type: expect.stringContaining("addAlgorithmNode"),
      }),
    )
  })

  it("adds a new input node when clicking the add button", () => {
    render(
      <Provider store={store}>
        <AlgorithmTreeView />
      </Provider>,
    )

    const dataAddButton = screen.getAllByLabelText("add")[0]
    fireEvent.click(dataAddButton)

    expect(store.dispatch).toHaveBeenCalledWith(
      expect.objectContaining({
        type: expect.stringContaining("addInputNode"),
      }),
    )
  })
})
