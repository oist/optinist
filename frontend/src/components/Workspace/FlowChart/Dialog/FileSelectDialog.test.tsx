import React, { createContext } from "react"
import { Provider } from "react-redux"

import { default as configureStore } from "redux-mock-store"

import { describe, it, beforeEach, jest, expect } from "@jest/globals"
import "@testing-library/jest-dom"
import { Store, AnyAction } from "@reduxjs/toolkit"
import { render, screen } from "@testing-library/react"

import { FileTreeItemLabel } from "components/Workspace/FlowChart/Dialog/FileSelectDialog"
import { AppDispatch } from "store/store"

// Create a mock context
type FileTreeActionsContextType = {
  onOpenDeleteDialog: (filePath: string, fileName: string) => void
}
const MockFileTreeActionsContext =
  createContext<FileTreeActionsContextType | null>(null)

const mockOnOpenDeleteDialog = jest.fn()

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
    mockOnOpenDeleteDialog.mockClear()
  })

  it("should render FileTreeItemLabel component", () => {
    render(
      <Provider store={store}>
        <MockFileTreeActionsContext.Provider
          value={{ onOpenDeleteDialog: mockOnOpenDeleteDialog }}
        >
          <FileTreeItemLabel
            multiSelect={true}
            fileType="all"
            shape={[100, 100]}
            label="testFile"
            isDir={false}
            checkboxProps={{ checked: false, onChange: jest.fn() }}
            filePath="testFile"
          />
        </MockFileTreeActionsContext.Provider>
      </Provider>,
    )

    // Check that the delete button exists
    const deleteButton = screen.getByTestId("DeleteIconBtn")
    expect(deleteButton).toBeTruthy()
  })

  it("should disable delete button if the file checkbox is checked", () => {
    render(
      <Provider store={store}>
        <MockFileTreeActionsContext.Provider
          value={{ onOpenDeleteDialog: mockOnOpenDeleteDialog }}
        >
          <FileTreeItemLabel
            multiSelect={true}
            fileType="all"
            shape={[100, 100]}
            label="testFile"
            isDir={false}
            checkboxProps={{ checked: true, onChange: jest.fn() }}
            filePath="testFile"
          />
        </MockFileTreeActionsContext.Provider>
      </Provider>,
    )

    const deleteButton = screen.getByTestId("DeleteIconBtn")
    expect(deleteButton.hasAttribute("disabled")).toBe(true)
  })

  it("should enable delete button if the file checkbox is not checked", () => {
    render(
      <Provider store={store}>
        <MockFileTreeActionsContext.Provider
          value={{ onOpenDeleteDialog: mockOnOpenDeleteDialog }}
        >
          <FileTreeItemLabel
            multiSelect={true}
            fileType="all"
            shape={[100, 100]}
            label="testFile"
            isDir={false}
            checkboxProps={{ checked: false, onChange: jest.fn() }}
            filePath="testFile"
          />
        </MockFileTreeActionsContext.Provider>
      </Provider>,
    )

    const deleteButton = screen.getByTestId("DeleteIconBtn")
    expect(deleteButton.hasAttribute("disabled")).toBe(false)
  })
})
