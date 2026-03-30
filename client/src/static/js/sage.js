const API_URL = process.env.API_URL

var TraceViewElement = require('./traceView.js').default;

window.customElements.define('trace-view', TraceViewElement);

var traceView = document.getElementById('traceView');

$('#mainTab a').on('click', function(e) {
  e.preventDefault()
  $(this).tab('show')
})

const resultLink = document.getElementById('link-results')

const submitButton = document.getElementById('btn-submit')
submitButton.addEventListener('click', showUpload)
const exampleButton = document.getElementById('btn-example')
exampleButton.addEventListener('click', showExample)

const inputFile = document.getElementById('inputFile')
const targetFastaFile = document.getElementById('targetFileFasta')
const targetChromatogramFile = document.getElementById('targetFileChromatogram')
const targetGenomes = document.getElementById('target-genome')
const targetTabs = document.getElementById('target-tabs')
const resultInfo = document.getElementById('result-info')
const resultError = document.getElementById('result-error')
const TRIM_MIN = 1
const TRIM_MAX = 1000

const GENOME_VALUE_RE = /^[A-Za-z0-9._-]{1,128}$/

function showValidationError(msg) {
  showElement(resultError)
  resultError.querySelector('#error-message').textContent = msg
}

function showExample() {
  run("example")
}

function showUpload() {
  run("data")
}

function run(stat) {
  resultLink.click()
  hideElement(resultError)
  traceView.deleteContent()

  const formData = new FormData()
  const lTrimRaw = document.getElementById('leftTrim').value
  const rTrimRaw = document.getElementById('rightTrim').value
  const lTrim = Number.parseInt(lTrimRaw, 10)
  const rTrim = Number.parseInt(rTrimRaw, 10)
  if ( (!Number.isInteger(lTrim)) || (lTrim < TRIM_MIN) || (lTrim > TRIM_MAX) ) {
    showValidationError(`Left trim must be an integer between ${TRIM_MIN} and ${TRIM_MAX}.`)
    return
  }
  if ( (!Number.isInteger(rTrim)) || (rTrim < TRIM_MIN) || (rTrim > TRIM_MAX) ) {
    showValidationError(`Right trim must be an integer between ${TRIM_MIN} and ${TRIM_MAX}.`)
    return
  }    
  formData.append('leftTrim', lTrim)
  formData.append('rightTrim', rTrim)
  if (stat === "example") {
    formData.append('showExample', 'showExample')
  } else {
    if ( (!inputFile.files) || (inputFile.files.length === 0) ) {
      showValidationError('Please select a chromatogram file.')
      return
    }
    formData.append('queryFile', inputFile.files[0])
    const target = targetTabs.querySelector('a.active').id
    if (target.startsWith('target-genome')) {
      const genomeVal = targetGenomes.querySelector('option:checked').value
      if (!GENOME_VALUE_RE.test(genomeVal)) {
        showValidationError('Invalid genome selection.')
	return
      }
      formData.append('genome', genomeVal)
    } else if (target.startsWith('target-fasta')) {
      if ( (!targetFastaFile.files) || (targetFastaFile.files.length === 0) ) {
        showValidationError('Please select a FASTA reference file.')
        return
      }
      formData.append('fastaFile', targetFastaFile.files[0])
    } else if (target.startsWith('target-chromatogram')) {
      if ( (!targetChromatogramFile.files) || (targetChromatogramFile.files.length === 0) ) {
	showValidationError('Please select a wildtype chromatogram file.')
	return
      }
      formData.append('chromatogramFile', targetChromatogramFile.files[0])
    }
  }
  
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
      traceView.deleteContent()
      showElement(resultError)
      resultError.querySelector('#error-message').textContent = errorMessage
    })
}

function handleSuccess(res) {
    hideElement(resultInfo)
    hideElement(resultError)
    traceView.displayData(res.data)
}

function showElement(element) {
  element.classList.remove('d-none')
}

function hideElement(element) {
  element.classList.add('d-none')
}



