import { Display } from './display.js';
import { isValid } from './utils.js';

export class StimulusDisplay extends Display {
    constructor(block, stims, feedbackPics, buttons) {
        super();
        this.clearContainer();
        this.insertBar();
        this.setTitle(`What's the orientation?<br> (Part ${block}/11)`);
       

        this.insertMain();
        
        // this._main.style.backgroundColor = 'blue';
        this._main.style.justifyContent = 'space-around';
        
        this.main_ITI = null;
        this.main_stimulus = null;

        this.pic_stimulus = null;

        this.quit_button = null;
        this.quit = null;

        this.stims = stims;
        this.stimsfiltered = stims.filter(isValid)
        this.feedbackPics = feedbackPics;

        this.feedback = null;

        this._create();
        this.setFoot(`Press <kbd>F</kbd> for <b>${buttons[0]}</b>, press <kbd>J</kbd> for <b>${buttons[1]}</b>`);
    }

    _create() {
        this._createMain_ITI(); 
        this._createMain_stimulus(); 
        this._createQuit();
        this._main.append(this.main_ITI, 
                        this.main_stimulus, 
            );
        this._container.insertBefore(this.quit, this._footer);
    }

    _createQuit() {
        if (this.quit === null){
            this.quit = document.createElement("div");
            this.quit.className = 'content-horizontal-right';
            this.quit.style.display = 'flex';
            this.quit.style.height = '10%';

            this.quit_button = document.createElement("button");
            this.quit_button.id = 'quit-button';
            this.quit_button.innerText = 'Quit';
            this.quit_button.style.display = 'inline-block';
            this.quit_button.style.textAlign = 'center';
            this.quit_button.style.fontSize = '1.8vh';
            this.quit_button.style.margin = "0px 20px 0px 0px";
            this.quit_button.style.borderColor = 'grey';
            this.quit_button.style.color = 'grey';
            this.quit_button.style.alignSelf = 'center';

            this.quit.append(this.quit_button);
        }
    }

    _createMain_ITI() {
        if (this.main_ITI === null) {
            this.main_ITI = document.createElement("div");
            this.main_ITI.id = 'main-iti';
            this.main_ITI.className = 'content-horizontal-center content-vertical-center';
            this.main_ITI.style.display = 'none';
            this.main_ITI.style.height = '50%';
            this.main_ITI.style.width = '30%';

            const fixation = document.createElement("p");
            fixation.id = "fixation";
            fixation.innerHTML = "";
            fixation.style.fontSize = "12vh";

            this.main_ITI.append(fixation);  
        }
    }

    _createMain_stimulus() {
        if (this.main_stimulus === null) {
            this.main_stimulus = document.createElement("div");
            this.main_stimulus.id = 'main-stimulus';
            this.main_stimulus.className = 'content-horizontal-center content-vertical-space-between';
            this.main_stimulus.style.display = 'flex';
            this.main_stimulus.style.order = 1;
            this.main_stimulus.style.margin = '10px';
            this.main_stimulus.style.height = '100%';
            this.main_stimulus.style.width = '50%';
            this.main_stimulus.style.flexDirection = 'column';
            // this.main_stimulus.style.backgroundColor = 'green';
            
            this.pic_stimulus = document.createElement("img");
            this.pic_stimulus.id = 'pic-stimulus';
            this.pic_stimulus.src = this.stimsfiltered[0];
            // this.pic_stimulus.style.height = '100%';
            this.pic_stimulus.style.margin = '5%';
            this.pic_stimulus.style.order = 1;
            this.pic_stimulus.style.visibility = 'visible';

            this.feedback = document.createElement("img");
            this.feedback.id = 'feedback';
            this.feedback.style.maxHeight = '25%';
            this.feedback.src = this.feedbackPics[0];
            this.feedback.style.order = 3;
            this.feedback.style.visibility = 'hidden';

            this.main_stimulus.append(this.pic_stimulus, this.feedback); 
            this._container.append(this.main_stimulus);
        }
    }


    showFixation(pic_source='', mirror) {
        this.pic_stimulus.src = pic_source;
        if (mirror) {
            this.pic_stimulus.style.transform = "scaleX(-1)";
        } else {
            this.pic_stimulus.style.transform = "scaleX(1)";
        }
        
        this.feedback.src = this.feedbackPics[0];

        this.feedback.style.visibility = 'visible';
        this.main_ITI.style.display = 'flex';
        this.main_stimulus.style.display = 'none';
    }

    showStimulus(val=50) {
        this.setBar(val);
        this.main_ITI.style.display = 'none';
        this.main_stimulus.style.display = 'flex';
    }

    showFeedback(isCorrect) {
        if (isCorrect) {
            this.feedback.src = this.feedbackPics[2];
        } else {
            this.feedback.src = this.feedbackPics[1];
        }

        this.feedback.style.visibility = 'visible';
    }

}



export function timeline(display, trialsmat, pics2nodes, studyID, what2press) {
    return new Promise( (resolve) => {
        let current_trial = -1;
        let isfinished = false;
        let sumCorrect = 0;
        let nResponses = 0;
        let results = [];
        nextTrial();

        function nextTrial() {
            if (isfinished) { 
                // display.clearContainer();
                display.clearMain();
                // display.clearFooter();

                let header = Object.keys(results[0]);
                results.unshift(header);

                let percentCorrect = Math.round(sumCorrect/nResponses * 100);

                display.setInstructions(`<p style="font-size:4vmin">Correct choices: <b>${percentCorrect}%</b></p>`);
                setTimeout( ()=>{resolve(results)}, 1500);
                // resolve(trials);     
            } else {
                runTrial();
            }
        }

        function runTrial() {

            current_trial++;
            let sourcenode = trialsmat[current_trial][0];
            let pic_source = pics2nodes[sourcenode];
            let mirror = trialsmat[current_trial][1];
            let distance = trialsmat[current_trial][2];
            let clustermember = trialsmat[current_trial][3];
            let keypress;
            let choice;
            let RT;
            let isCorrect;
            let timeout;

            display.showFixation(
                display.stims[pic_source], 
                mirror,
                );

            let trialonset = null;
            let trialoffset = null;

            setTimeout(stim_on, 500);
        
            function stim_on() {
                
                const progress = ((current_trial+1)/trialsmat.length)*100;
                display.showStimulus(progress);
                trialonset = new Date().getTime();
                document.addEventListener('keydown', feedback, {once: true});
                timeout = setTimeout(offset, 2000);
                isCorrect = 0;
                choice = 0;
                RT = 'NA';
                keypress = 'NA';
            }

            function feedback(event) {
                
                clearTimeout(timeout);
                trialoffset = new Date().getTime();
                RT = trialoffset - trialonset;
                keypress = event.code;
                switch (keypress) {
                    case 'KeyF':
                        choice = 1;
                        nResponses = nResponses + 1;
                        if (0 === mirror & what2press[0] === 'normal') {
                            isCorrect = 1;
                        } else if (1 === mirror & what2press[0] === 'mirrored') {
                            isCorrect = 1;
                        }
                        display.showFeedback(isCorrect);
                        break;
                    case 'KeyJ':
                        choice = 2;
                        nResponses = nResponses + 1;
                        if (0 === mirror & what2press[1] === 'normal') {
                            isCorrect = 1;
                        } else if (1 === mirror & what2press[1] === 'mirrored') {
                            isCorrect = 1;
                        }
                        display.showFeedback(isCorrect);
                        break;
                }
                setTimeout(offset, 2000 - RT);
                
            }

            function offset() {
                document.removeEventListener('keydown', feedback, {once: true});
                sumCorrect += isCorrect;
                

                let result = new Result(current_trial+1,
                                        sourcenode,
                                        pic_source,
                                        mirror,
                                        keypress, 
                                        choice,
                                        RT,
                                        isCorrect,
                                        distance, 
                                        clustermember,
                                        studyID);
                results.push(result);
                nextTrial();
            }

            // if this was the last trial
            if (current_trial === trialsmat.length-1){
                isfinished = true;
            }
        }
        class Result {
            constructor(trialid, 
                        sourcenode, sourcepic, 
                        mirror,
                        keypress,
                        choice,
                        RT,
                        isCorrect, 
                        distance,
                        clustermember,
                        graph = studyID,
                        taskid = 'SEQ',) {
                this.trialid = trialid;
                this.taskid = taskid;
                this.graph = graph;
                this.distance = distance;
                this.clustermember = clustermember;
                this.sourcenode = sourcenode;
                this.sourcepic = sourcepic;
                this.mirror = mirror;
                this.keypress = keypress;
                this.choice = choice; // left or right
                this.RT = RT; // reaction time
                this.isCorrect = isCorrect; 
            }
        }
    })
}