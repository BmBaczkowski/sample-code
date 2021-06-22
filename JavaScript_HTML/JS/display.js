export class Display {
    constructor(className) {
        this._container = document.getElementById("container");
        this._bar = document.getElementById("bar");
        this._header = document.getElementById("header");
        this._main = document.getElementById("main");
        this._footer = document.getElementById("footer");
        this.quit_button = document.getElementById("quit-button");
        this.start_button = document.getElementById("start-button");
        this.clearContainer();
        this.insertMain(className);
    }

    insertBar() {
        this._bar = document.createElement("div");
        this._bar.id = 'bar';
        this._bar.className = 'progress';
        this._bar.style.backgroundColor = 'green';
        this._bar.style.width = '0%';
        if (this._header === null) {this.insertHeader()};
        this._container.insertBefore(this._bar, this._header);
    }

    insertMain(className="content-horizontal-center content-vertical-center") {
        if (this._main === null) {
            this._main = document.createElement("main");
            this._main.id = "main";
            this._main.className = className;
            this._main.style.minHeight = '65%';
            this._main.style.flexDirection = 'column';
            this._container.appendChild(this._main);

            console.log("Main created");
        } else {
            console.log("Main already exists");
        }
    }

    setInstructions(str) {
        if (this._main === null) {
            this.insertMain();
        }
        this._main.innerHTML= str;
        this._main.style.fontSize = '2.5vh';
    }

    insertButtons(txtleft='Start', txtright='Quit') {
        if (this._footer === null) {
            this.insertFooter('content-vertical-center content-horizontal-space-between');
        } 
        this.insertFooter('content-vertical-center content-horizontal-space-between');

        this.placeholder = document.createElement("button");
        this.placeholder.id = 'start-button';
        this.placeholder.textContent = txtleft;
        this.placeholder.style.display = 'inline-block';
        this.placeholder.style.fontSize = '1.8vh';
        this.placeholder.style.textAlign = 'center';
        this.placeholder.style.padding = "1vmin 3vmin";
        this.placeholder.style.margin = "2vmin";
        this.placeholder.style.borderColor = 'green';
        this.placeholder.style.color = 'green';
        this.placeholder.style.visibility = 'hidden';

        this.start_button = document.createElement("button");
        this.start_button.id = 'start-button';
        this.start_button.textContent = txtleft;
        this.start_button.style.display = 'inline-block';
        this.start_button.style.fontSize = '1.8vh';
        this.start_button.style.textAlign = 'center';
        this.start_button.style.padding = "1vmin 3vmin";
        this.start_button.style.margin = "2vmin";
        this.start_button.style.borderColor = 'green';
        this.start_button.style.color = 'green';
        
        this.quit_button = document.createElement("button");
        this.quit_button.id = 'quit-button';
        this.quit_button.innerText = txtright;
        this.quit_button.style.display = 'inline-block';
        this.quit_button.style.textAlign = 'center';
        this.quit_button.style.padding = "1vmin 3vmin";
        this.quit_button.style.fontSize = '1.8vh';
        this.quit_button.style.margin = "2vmin";  
        // this.quit_button.style.borderColor = 'red';
        // this.quit_button.style.color = 'red';
        
        this._footer.append(this.placeholder, this.start_button, this.quit_button);
    }

    insertHeader(className="content-horizontal-center content-vertical-center") {
        if (this._header === null) {
            this._header = document.createElement("header");
            this._header.id = "header";
            this._header.className = className;
            this._header.style.maxHeight = '20%';
            // if (this._main === null) {this.insertMain()};
            this._container.insertBefore(this._header, this._main);
            console.log("Header created");
        } else {
            console.log("Header already exists");
        }
    }

    insertFooter(className="content-horizontal-center content-vertical-center") {
        if (this._footer === null) {
            this._footer = document.createElement("footer");
            this._footer.id = "footer";
            this._footer.className = className;
            this._footer.style.maxHeight = '15%';
            // if (this._main === null) {this.insertMain()};
            this._container.append(this._footer);
            console.log("Footer created");  
        } else {
            console.log("Footer already exists");
        }
    }

    showFixation() {
        if (this._main === null) {
            this.insertMain();
        }
        
        const fixation = document.createElement("p");
        fixation.innerHTML = '+'
        fixation.id = "fixation";
        fixation.innerHTML = "+";
        fixation.style.fontSize = "12vh";

        this._main.append(fixation);  
    }

    setBar(val) {
        if (this._bar === null) {
            this.insertBar();
            this._bar.style.width = `${val}%`;
        } else {
            this._bar.style.width = `${val}%`;
        }
    }

    setTitle(str="") {
        if (this._header === null) {
            this.insertHeader('content-horizontal-center');
        } else {
            this._header.textContent = '';
        }
        
        let title = document.createElement("h3");
        title.style.fontSize = '3vmin';
        title.innerHTML = str;
        this._header.appendChild(title);
    }

    setFoot(str="") {
        if (this._footer === null) {
            this.insertFooter();
        } else {
            this._footer.textContent = '';
        }
        let foot = document.createElement("p");
        foot.innerHTML = str;
        foot.style.maxWidth = '80%';
        foot.style.fontSize = '3vmin';
        this._footer.appendChild(foot);
    }

    clearMain() {
        if (this._main !== null) {
            this._main.textContent = '';
        }
    }

    clearHeader() {
        if (this._header !== null) {
            this._header.textContent = '';
        }

    }

    clearFooter() {
        if (this._footer !== null) {
            this._footer.textContent = '';
        }
    }

    clearContainer() {
        this._container.textContent = '';
        this._header = null;
        this._main = null;
        this._footer = null;
    }
}