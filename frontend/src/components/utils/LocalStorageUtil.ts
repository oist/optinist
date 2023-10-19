import React from "react"

export function useLocalStorage<T>(
  storageKey: string,
  initialValue: T,
  parseFn: (value: any) => T,
): [T, Dispatch<SetStateAction<T>>] {
  const [value, setValue] = useState(() => {
    const savedStr = localStorage.getItem(storageKey)
    if (savedStr != null) {
      return parseFn(JSON.parse(savedStr))
    }
    return initialValue
  })
  useEffect(() => {
    localStorage.setItem(storageKey, JSON.stringify(value))
  }, [value, storageKey])
  return [value, setValue]
}
