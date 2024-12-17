import React, { useState } from "react"
import { Provider } from "react-redux"

import { SnackbarProvider } from "notistack"
import configureStore from "redux-mock-store"

import { describe, it, beforeEach } from "@jest/globals"
import { render, screen } from "@testing-library/react"
import { userEvent } from "@testing-library/user-event"

import { RunButtons } from "components/Workspace/FlowChart/Buttons/RunButtons"
import { RUN_BTN_OPTIONS } from "store/slice/Pipeline/PipelineType"

const mockStore = configureStore([])

describe("RunButtons component", () => {
  let store: ReturnType<typeof mockStore>

  const mockProps = {
    uid: "test-uid",
    runDisabled: false,
    filePathIsUndefined: false,
    algorithmNodeNotExist: false,
    handleCancelPipeline: jest.fn(),
    handleRunPipeline: jest.fn(),
    handleRunPipelineByUid: jest.fn(),
  }

  beforeEach(() => {
    store = mockStore({
      pipeline: {
        run: {
          status: "StartUninitialized",
        },
      },
      currentPipeline: {
        uid: "test-uid",
      },
      runBtn: RUN_BTN_OPTIONS.RUN_NEW,
    })
  })

  it("renders correctly and triggers handleRunPipelineByUid on Run button click", async () => {
    render(
      <Provider store={store}>
        <SnackbarProvider>
          <RunButtons status={"StartUninitialized"} {...mockProps} />
        </SnackbarProvider>
      </Provider>,
    )

    // Find the button using the PlayArrow icon's testid
    const runButtonIcon = screen.getByTestId("PlayArrowIcon")
    // Get the parent button from the icon
    const runButton = runButtonIcon.closest("button")

    if (runButton) {
      await userEvent.click(runButton)

      // Check if handleRunPipelineByUid was called
      expect(mockProps.handleRunPipelineByUid).toHaveBeenCalledTimes(1)
    } else {
      throw new Error("Run button not found")
    }
  })

  it("disables the button after it is clicked", async () => {
    // Use state to simulate the button becoming disabled after a click
    const WrapperComponent = () => {
      const [runDisabled, setRunDisabled] = useState(false)

      const updatedProps = {
        ...mockProps,
        runDisabled,
        handleRunPipelineByUid: () => {
          // Simulate disabling the button after clicking
          setRunDisabled(true)
        },
      }

      return (
        <Provider store={store}>
          <SnackbarProvider>
            <RunButtons status={"StartUninitialized"} {...updatedProps} />
          </SnackbarProvider>
        </Provider>
      )
    }

    render(<WrapperComponent />)

    // Find the button using the PlayArrow icon's testid
    const runButtonIcon = screen.getByTestId("PlayArrowIcon")

    // Get the parent button from the icon
    const runButton = runButtonIcon.closest("button")

    if (runButton) {
      // Initially, the button should not be disabled
      expect(runButton).not.toBeDisabled()

      // Click the button
      await userEvent.click(runButton)

      // Verify that the button is disabled after the click
      expect(runButton).toBeDisabled()
    } else {
      throw new Error("Run button not found")
    }
  })
})
