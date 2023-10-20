import {
  useCallback,
  DependencyList,
  MouseEvent as ReactMouseEvent,
} from "react"

export function useMouseDragHandler(
  onMouseDown: (event: ReactMouseEvent) => {
    onMouseMove: (this: Document, event: MouseEvent) => void
    onMouseUp: (this: Document, event: MouseEvent) => void
  },
  dependencies: DependencyList,
) {
  return useCallback(
    (event: ReactMouseEvent) => {
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
