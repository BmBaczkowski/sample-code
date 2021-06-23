import { Display } from './display.js';
import { shuffle, getVect } from './utils.js';


export class StimulusDisplay extends Display {
    constructor(stims) {
        super();
        this.clearContainer();
        this.setTitle('Task');
        this.insertMain("content-horizontal-center content-vertical-space-around");
        this.insertButtons('Start', 'Quit');
        
        // this._main.style.backgroundColor = 'blue';
        this._main.style.flexDirection = 'column';

        this._create();

        this.stims = [];

        let ivect = getVect(stims.length);
        ivect.shift();
        ivect = shuffle(ivect);

        // let pics2nodes = pics2nodesObj[graph];
        for (let i of ivect) {
            this.stims[i] = document.createElement("img");
            this.stims[i].id = `img-${i}`;  
            this.stims[i].src = stims[i];
            this.stims[i].style.height = '50%';
            this.stims[i].style.width = 'auto';       

            this.stims[i].style.borderWidth = 'thin';
            this.stims[i].style.margin = '2px';
            this.stims[i].style.borderStyle= 'solid';
            this.stims[i].style.borderRadius = '10%';
            this.stims[i].style.borderColor = 'black';
            
            this.img_div.append(this.stims[i]);
            
        }

    }

    _create() {
        this._createInstDiv();
        this._createImgDiv();
    }

    _createInstDiv() {
        this.inst_div = document.createElement("div");
        this.inst_div.id = 'instructions';
        this.inst_div.className = 'content-horizontal-center';
        this.inst_div.style.width = '80%';
        this.inst_div.style.height = '50%';
        // this.inst_div.style.backgroundColor = 'lightgrey';
        this.inst_div.style.order = 1;

        this._main.append(this.inst_div);
    }

    _createImgDiv() {
        this.img_div = document.createElement("div");
        this.img_div.id = 'img-div';
        this.img_div.className = 'content-horizontal-space-around content-vertical-center';
        this.img_div.style.width = '60%';
        this.img_div.style.maxHeight = '40%';
        // this.img_div.style.backgroundColor = 'red';
        this.img_div.style.order = 2;

        this._main.append(this.img_div);
    }

    setInstr(instr = '') {
        this.inst_div.innerHTML = instr;
    }

}
