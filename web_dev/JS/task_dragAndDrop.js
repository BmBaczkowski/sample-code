import { Display } from './display.js';
import { shuffle, getVect } from './utils.js';


export class StimulusDisplay extends Display {
    constructor(stims) {
        super();
        this.clearContainer();
        this.setTitle('Arrange objects in the grey area');
        this.insertMain("content-horizontal-center content-vertical-space-around");
        
        // this._main.style.backgroundColor = 'blue';
        this._main.style.flexDirection = 'column';

        this.submit_button = null;
        this.reset_button = null;
        this.quit_button = null;


        this._create();

        this.stims = [];

        let ivect = getVect(stims.length);
        ivect.shift();
        ivect = shuffle(ivect);

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
        this._createArena(); 
        this._createImgDiv();
        this._createButtons(); 
    }

    _createImgDiv() {
        this.img_div = document.createElement("div");
        this.img_div.id = 'img-div';
        this.img_div.className = 'content-horizontal-space-around content-vertical-center';
        this.img_div.style.width = '60%';
        this.img_div.style.height = '20%';
        // this.img_div.style.backgroundColor = 'red';
        this.img_div.style.order = 2;

        this._main.append(this.img_div);
    }

    _createArena() {
        this.arena_div = document.createElement("div");
        this.arena_div.id = 'arena';
        this.arena_div.className = 'circle-container';
        this.arena_div.style.width = `${70+10}vmin`;
        this.arena_div.style.height = '70vmin';
        this.arena_div.style.backgroundColor = 'lightgrey';
        this.arena_div.style.order = 1;

        this._main.append(this.arena_div);
    }


    _createButtons() {
        this.insertFooter()

        this.submit_button = document.createElement("button");
        this.submit_button.id = 'submit-button';
        this.submit_button.textContent = 'Submit';
        this.submit_button.style.display = 'inline-block';
        this.submit_button.style.fontSize = '1.8vh';
        this.submit_button.style.textAlign = 'center';
        this.submit_button.style.padding = "1vmin 3vmin";
        this.submit_button.style.margin = "2vmin";


        this.reset_button = document.createElement("button");
        this.reset_button.id = 'reset-button';
        this.reset_button.textContent = 'Reset';
        this.reset_button.style.display = 'inline-block';
        this.reset_button.style.fontSize = '1.8vh';
        this.reset_button.style.textAlign = 'center';
        this.reset_button.style.padding = "1vmin 3vmin";
        this.reset_button.style.margin = "2vmin";  
        
        this.quit_button = document.createElement("button");
        this.quit_button.id = 'quit-button';
        this.quit_button.innerText = 'Quit';
        this.quit_button.style.display = 'inline-block';
        this.quit_button.style.textAlign = 'center';
        this.quit_button.style.padding = "1vmin 3vmin";
        this.quit_button.style.fontSize = '1.8vh';
        this.quit_button.style.margin = "2vmin";  
        // this.quit_button.style.borderColor = 'grey';
        // this.quit_button.style.color = 'grey';

        this._footer.append(this.submit_button, this.reset_button, this.quit_button);
    }

}
