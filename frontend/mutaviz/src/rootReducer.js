import { combineReducers } from 'redux'
import appReducer from './reducers/appReducer'
import formReducer from './reducers/formReducer'

export const rootReducer = combineReducers({
  main: appReducer,
  form: formReducer
})
