type StateAction =
  | { type: 'AlgoSelect'; value: string }
  | { type: 'ParamUpdate'; value: number | string; param: string }

export default StateAction
