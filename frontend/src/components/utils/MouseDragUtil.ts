import React from 'react'

export function useMouseDragHandler(
  onMouseDown: (event: React.MouseEvent) => {
    onMouseMove: (this: Document, event: MouseEvent) => void
    onMouseUp: (this: Document, event: MouseEvent) => void
  },
  dependencies: React.DependencyList,
) {
  return React.useCallback(
    (event: React.MouseEvent) => {
      const { onMouseMove, onMouseUp } = onMouseDown(event)
      document.addEventListener('mousemove', onMouseMove)
      document.addEventListener(
        'mouseup',
        (event) => {
          document.removeEventListener('mousemove', onMouseMove)
          onMouseUp.call(document, event)
        },
        { once: true },
      )
    },
    // eslint-disable-next-line react-hooks/exhaustive-deps
    [dependencies],
  )
}
