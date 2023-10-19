import React from "react"

export function useMouseDragHandler(
  onMouseDown: (event: MouseEvent) => {
    onMouseMove: (this: Document, event: MouseEvent) => void
    onMouseUp: (this: Document, event: MouseEvent) => void
  },
  dependencies: DependencyList,
) {
  return useCallback(
    (event: MouseEvent) => {
      const { onMouseMove, onMouseUp } = onMouseDown(event)
      document.addEventListener("mousemove", onMouseMove)
      document.addEventListener(
        "mouseup",
        (event) => {
          document.removeEventListener("mousemove", onMouseMove)
          onMouseUp.call(document, event)
        },
        { once: true },
      )
    },
    // eslint-disable-next-line react-hooks/exhaustive-deps
    [dependencies],
  )
}
