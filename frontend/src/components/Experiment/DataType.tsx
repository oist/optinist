export function createData(
  date: string,
  status: boolean,
  progress: string,
  name: string,
) {
  return {
    date,
    status,
    progress,
    name,
    details: [
      { function: 'suite2p_file_convert', success: true },
      { function: 'suite2p_registration', success: false },
    ],
  }
}
