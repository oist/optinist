import { describe, expect, test } from "@jest/globals"

import { splitValue } from "components/common/ParamFormItemCreator"

describe("splitValue", () => {
  test.each([
    ["0,1", [0, 1]],
    [" t , b ", ["t", "b"]],
    ["1,,2", [1, 2]],
    ["-1, 1e3", [-1, 1000]],
    ["0, tbc", [0, "tbc"]],
  ])("%p parses to %p", (input, expected) => {
    expect(splitValue(input)).toEqual(expected)
  })
})
