import React from "react"
import { Provider } from "react-redux"

import configureStore from "redux-mock-store"

import { Store, AnyAction } from "@reduxjs/toolkit"
import { render, screen, fireEvent } from "@testing-library/react"

import { TreeItemLabel } from "components/Workspace/FlowChart/Dialog/FileSelectDialog"
import { deleteFileTree } from "store/slice/FilesTree/FilesTreeAction"
import { AppDispatch } from "store/store"

jest.mock("store/slice/FilesTree/FilesTreeAction", () => ({
  deleteFileTree: jest.fn(),
}))

const mockStore = configureStore<
  Partial<{ workspace: { currentWorkspace: { workspaceId?: number } } }>,
  AppDispatch
>([])

describe("TreeItemLabel Component", () => {
  let store: Store<unknown, AnyAction>

  beforeEach(() => {
    store = mockStore({
      workspace: {
        currentWorkspace: {
          workspaceId: 123,
        },
      },
    })
    store.dispatch = jest.fn()
  })

  it("should dispatch deleteFileTree action when delete is confirmed", () => {
    render(
      <Provider store={store}>
        <TreeItemLabel
          fileType="image"
          shape={[100, 100]}
          label="testFile"
          isDir={false}
          checkboxProps={{ checked: false, onChange: jest.fn() }}
          setSelectedFilePath={jest.fn()}
          selectedFilePath={""}
        />
      </Provider>,
    )

    // Simulate opening the delete confirmation dialog
    const deleteIcon = screen.getByTestId("DeleteIcon")
    fireEvent.click(deleteIcon)

    // Confirm the delete action
    const confirmButton = screen.getByRole("button", { name: /yes/i })
    fireEvent.click(confirmButton)

    expect(store.dispatch).toHaveBeenCalledWith(
      deleteFileTree({
        workspaceId: 123,
        fileName: "testFile",
        fileType: "image",
      }),
    )
  })

  it("should not dispatch deleteFileTree action if confirmation clicked no", () => {
    render(
      <Provider store={store}>
        <TreeItemLabel
          fileType="image"
          shape={[100, 100]}
          label="testFile"
          isDir={false}
          checkboxProps={{ checked: false, onChange: jest.fn() }}
          setSelectedFilePath={jest.fn()}
          selectedFilePath={""}
        />
      </Provider>,
    )

    // Simulate opening the delete confirmation dialog
    const deleteIcon = screen.getByTestId("DeleteIcon")
    fireEvent.click(deleteIcon)

    // Confirm the delete action
    const confirmButton = screen.getByRole("button", { name: /no/i })
    fireEvent.click(confirmButton)

    expect(store.dispatch).not.toHaveBeenCalledWith(deleteFileTree)
  })
})
