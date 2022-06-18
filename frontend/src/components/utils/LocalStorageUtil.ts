import React from 'react'

export function useLocalStorage<T>(
  storageKey: string,
  initialValue: T,
  parseFn: (value: any) => T,
): [T, React.Dispatch<React.SetStateAction<T>>] {
  const [value, setValue] = React.useState(() => {
    const savedStr = localStorage.getItem(storageKey)
    if (savedStr != null) {
      return parseFn(JSON.parse(savedStr))
    }
    return initialValue
  })
  React.useEffect(() => {
    localStorage.setItem(storageKey, JSON.stringify(value))
  }, [value, storageKey])
  return [value, setValue]
}
