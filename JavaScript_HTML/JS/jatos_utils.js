export function abortBox(id) {

    let response = prompt(`To quit the study, click OK and tell us whether\nwe can keep your data collected in this session.\nType "Yes" or "No".`, 'Yes');

    if (response !== null) {
        if (typeof jatos.batchSession.find(`/${id}/studyID`) === 'undefined') {
            jatos.batchSession.add(`/${id}`, {studyID: 'xxx'})
            .then( () => {
                if (response.toLowerCase().includes('yes') === true) {
                    jatos.endStudyAjax(false, `Prolific id: ${id}. Participant aborted by pressing abort button, but OK to keep the data.`, endpage);
                } else {
                    jatos.abortStudyAjax(`Prolific id: ${id}. Participant aborted by pressing abort button.`, endpage);
                }

            })
        } else {
            jatos.batchSession.replace(`/${id}/studyID`, 'xxx')
            .then( () => {
                if (response.toLowerCase().includes('yes') === true) {
                    jatos.endStudyAjax(false, `Prolific id: ${id}. Participant aborted by pressing abort button, but OK to keep the data.`, endpage);
                } else {
                    jatos.abortStudyAjax(`Prolific id: ${id}. Participant aborted by pressing abort button.`, endpage);
                }

            })
        }
    }

    function endpage() {
        document.body.innerHTML = `
        <div style="text-align:center">
        <p>
            The session was cancelled.
        </p>
        </div>
        `;
    }
}
