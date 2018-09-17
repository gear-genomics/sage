const API_URL = process.env.API_URL

var traceView = require('./traceView');

$('#mainTab a').on('click', function(e) {
  e.preventDefault()
  $(this).tab('show')
})

const resultLink = document.getElementById('link-results')

const submitButton = document.getElementById('btn-submit')
submitButton.addEventListener('click', function() {
  resultLink.click()
  run()
})

const exampleButton = document.getElementById('btn-example')
exampleButton.addEventListener('click', showExample)

const inputFile = document.getElementById('inputFile')
const targetFastaFile = document.getElementById('targetFileFasta')
const targetChromatogramFile = document.getElementById('targetFileChromatogram')
const targetGenomes = document.getElementById('target-genome')
const targetTabs = document.getElementById('target-tabs')
const resultContainer = document.getElementById('result-container')
const resultInfo = document.getElementById('result-info')
const resultError = document.getElementById('result-error')

// TODO client-side validation
function run() {
  const formData = new FormData()
  formData.append('queryFile', inputFile.files[0])
  const target = targetTabs.querySelector('a.active').id

  if (target.startsWith('target-genome')) {
    const genome = targetGenomes.querySelector('option:checked').value
    formData.append('genome', genome)
  } else if (target.startsWith('target-fasta')) {
    formData.append('fastaFile', targetFastaFile.files[0])
  } else if (target.startsWith('target-chromatogram')) {
    formData.append('chromatogramFile', targetChromatogramFile.files[0])
  }
  
  traceView.deleteContent()
  hideElement(resultContainer)
  hideElement(resultError)
  showElement(resultInfo)

  axios
    .post(`${API_URL}/upload`, formData)
    .then(res => {
	if (res.status === 200) {
          handleSuccess(res.data)
      }
    })
    .catch(err => {
      let errorMessage = err
      if (err.response) {
        errorMessage = err.response.data.errors
          .map(error => error.title)
          .join('; ')
      }
      hideElement(resultInfo)
      showElement(resultError)
      resultError.querySelector('#error-message').textContent = errorMessage
    })
}

async function handleSuccess(res) {
    hideElement(resultInfo)
    hideElement(resultError)
    showElement(resultContainer)
    traceView.displayData(res.data)
}

function showExample() {
    resultLink.click()
    const formData = new FormData()
    formData.append('showExample', 'showExample')
    hideElement(resultContainer)
    hideElement(resultError)
    showElement(resultInfo)
    axios
	.post(`${API_URL}/upload`, formData)
	.then(res => {
	    if (res.status === 200) {
		handleSuccess(res.data)
	    }
	})
	.catch(err => {
	    let errorMessage = err
	    if (err.response) {
		errorMessage = err.response.data.errors
		    .map(error => error.title)
		    .join('; ')
	    }
	    hideElement(resultInfo)
	    showElement(resultError)
	    resultError.querySelector('#error-message').textContent = errorMessage
	})
}

function showElement(element) {
  element.classList.remove('d-none')
}

function hideElement(element) {
  element.classList.add('d-none')
}



