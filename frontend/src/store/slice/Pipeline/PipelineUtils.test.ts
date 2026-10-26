import { OptionsObject, SnackbarKey } from "notistack"

import { describe, expect, jest, test } from "@jest/globals"

import { handleWorkflowYamlError } from "store/slice/Pipeline/PipelineUtils"

const enqueue = () =>
  jest.fn((_message: string, _options?: OptionsObject): SnackbarKey => "key")

describe("handleWorkflowYamlError", () => {
  test("a 422 shows the backend detail with the FAQ action", () => {
    const snackbar = enqueue()
    const detail =
      "Workflow yaml error, see FAQ: unknown parameter 'transpose' for dpca; delete the node, add it again and rerun"
    handleWorkflowYamlError(
      { response: { status: 422, data: { detail } } },
      snackbar,
    )
    expect(snackbar.mock.calls[0][0]).toBe(`${detail}\n`)
    expect(snackbar.mock.calls[0][1]).toHaveProperty("action")
  })

  test("a 422 without a detail falls back to the generic text", () => {
    const snackbar = enqueue()
    handleWorkflowYamlError({ response: { status: 422 } }, snackbar)
    expect(snackbar.mock.calls[0][0]).toBe("Workflow yaml error, see FAQ\n")
  })

  test("any other status is the generic failure", () => {
    const snackbar = enqueue()
    handleWorkflowYamlError({ response: { status: 500 } }, snackbar)
    expect(snackbar.mock.calls[0][0]).toBe("Failed to Run workflow")
  })
})
