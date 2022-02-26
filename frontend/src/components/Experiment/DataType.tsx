export function createData(
  date: string,
  name: string,
  status: boolean,
  progress: number,
) {
  return {
    date,
    name,
    status,
    progress,
    details: [
      { function: 'suite2p_file_convert', success: true },
      { function: 'suite2p_registration', success: false },
    ],
  }
}
