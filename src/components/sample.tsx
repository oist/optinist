import React from 'react'
import ReactDOM from 'react-dom'
import FlexLayout, { IJsonModel, Model } from 'flexlayout-react'

const json: IJsonModel = {
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
            name: 'FX',
            component: 'button',
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
            name: 'FI',
            component: 'button',
          },
        ],
      },
    ],
  },
}

interface MyState {
  model: Model
}

class Sample extends React.Component<{}, MyState> {
  constructor(props: object) {
    super(props)
    this.state = { model: FlexLayout.Model.fromJson(json) }
  }

  factory = (node: any) => {
    var component = node.getComponent()
    if (component === 'button') {
      return <button>{node.getName()}</button>
    } else {
      return null
    }
  }

  render() {
    return <FlexLayout.Layout model={this.state.model} factory={this.factory} />
  }
}

ReactDOM.render(<Sample />, document.getElementById('container'))
