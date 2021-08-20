import './App.css'
import { Layout, Model, TabNode, IJsonModel } from 'flexlayout-react'
import 'flexlayout-react/style/light.css'
import BasicFlow from 'components/flow'

var json: IJsonModel = {
  global: {},
  layout: {
    type: 'row',
    weight: 100,
    children: [
      {
        type: 'tabset',
        weight: 50,
        selected: 0,
        children: [
          {
            type: 'tab',
            name: 'One',
            component: 'flow',
          },
        ],
      },
      {
        type: 'tabset',
        weight: 50,

        selected: 0,
        children: [
          {
            type: 'tab',
            name: 'Two',

            component: 'grid',
          },
        ],
      },
    ],
  },
  borders: [
    {
      type: 'border',
      location: 'bottom',
      size: 100,
      children: [
        {
          type: 'tab',
          name: 'four',
          component: 'grid',
        },
      ],
    },
    {
      type: 'border',
      location: 'left',
      size: 100,
      children: [],
    },
  ],
}

const model = Model.fromJson(json)

function App() {
  const factory = (node: TabNode) => {
    var component = node.getComponent()
    if (component === 'button') {
      return <button>{node.getName()}</button>
    } else if (component == 'flow') {
      return <BasicFlow />
    } else {
      return null
    }
  }

  return <Layout model={model} factory={factory} />
}

export default App
