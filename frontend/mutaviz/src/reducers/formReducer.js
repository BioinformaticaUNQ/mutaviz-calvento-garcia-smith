const initialState = {
  fasta: "",
  sequence: "",
  fastaValid: false,
  goToLine: null,
  mutations: [],
  proteineStartPosition: null,
  selectedAminoacid: {
    position: null,
    from: ""
  },
}

const reducer = (state = initialState, action) => {
  return state
}
export default reducer