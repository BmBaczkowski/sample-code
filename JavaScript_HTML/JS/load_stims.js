const fetch_retry = async (url, n) => {

    let response = await fetch(url);
    if (response.ok === true) {
        return response;
    } else {
        if (n === 1) {return response};
        return await fetch_retry(url, n - 1);
    } 
};

export async function catchStims(picnumbers) {
    let stims = [];
    let str = '';
    let response = null;
    let blob = null;

    // remove null
    picnumbers = picnumbers.slice(1);

    for (let i of picnumbers) {
        str = `./stimuli/${('00'+i).slice(-2)}.png`;

        // response = await fetch(str);
        response = await fetch_retry(str, 5);
        if (response.ok === false) {
            return false;
        }
        blob = await response.blob();
        stims[i] = URL.createObjectURL(blob);
    }
    return stims;
}

export async function catchFeedback() {
    let feedback = [];
    let str = '';
    let response = null;
    let blob = null;

    str = `./stimuli/feedback/placeholder.png`;
    response = await fetch_retry(str, 5);
    if (response.ok === false) {
        return false;
    }
    blob = await response.blob();
    feedback[0] = URL.createObjectURL(blob);

    str = `./stimuli/feedback/incorrect.png`;
    response = await fetch_retry(str, 5);
    if (response.ok === false) {
        return false;
    }
    blob = await response.blob();
    feedback[1] = URL.createObjectURL(blob);

    str = `./stimuli/feedback/correct.png`;
    response = await fetch_retry(str, 5);
    if (response.ok === false) {
        return false;
    }
    blob = await response.blob();
    feedback[2] = URL.createObjectURL(blob);
    return feedback;
}

export async function catchFrame() {
    let str = `./attachments/associative_frame.png`;
    let response = await fetch_retry(str, 5);
    if (response.ok === false) {
        return false;
    }
    let blob = await response.blob();
    let frame = URL.createObjectURL(blob);
    return frame;
}