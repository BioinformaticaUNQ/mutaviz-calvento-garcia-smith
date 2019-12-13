import React from 'react'
import { connect } from 'react-redux'
import { LOADING, FORM, RESULT } from '../../reducers/appReducer'
import Loading from '../Loading'
import ConnectedForm from '../Form'
import Result from '../Result'

const contentFromPageState = (state) => {
  switch (state) {
    case LOADING: return <Loading />
    case FORM: return <ConnectedForm />
    case RESULT: return <Result />
  }
}

const Body = ({ pageState }) => {
  return (
    <body>
      {contentFromPageState(pageState)}
    </body>
  )
}

const mapStateToProps = (state) => {
  return {
    pageState: state.main.pageState
  }
}

export default connect(mapStateToProps, null)(Body)