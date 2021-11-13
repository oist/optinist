export function arrayEqualityFn<T>(a: T[], b: T[]) {
  return a === b || (a.length === b.length && a.every((v, i) => v === b[i]))
}

export function twoDimarrayEqualityFn<T>(a: T[][], b: T[][]) {
  return (
    a === b ||
    (a.length === b.length && a.every((v, i) => arrayEqualityFn(v, b[i])))
  )
}
