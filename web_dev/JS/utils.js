export function getDate() {
    let d = new Date();
    let hour = '' + d.getHours();
    let minutes = '' + d.getMinutes();
    let month = '' + (d.getMonth() + 1);
    let day = '' + d.getDate();
    let year = d.getFullYear();

    if (hour.length < 2) 
        hour = '0' + hour;  
    if (minutes.length < 2) 
        minutes = '0' + minutes; 
    if (month.length < 2) 
        month = '0' + month;
    if (day.length < 2) 
        day = '0' + day;

    let time = [hour, minutes].join(':'); 
    let date = [day, month, year].join('-'); 
    let timedate = time + '//' + date;
    return timedate;
};

export function getVisibilityAPI() {
    // based on MDN
    // Set the name of the hidden property and the change event for visibility
    let hidden, visibilityChange; 
    if (typeof document.hidden !== "undefined") { // Opera 12.10 and Firefox 18 and later support 
    hidden = "hidden";
    visibilityChange = "visibilitychange";
    } else if (typeof document.msHidden !== "undefined") {
    hidden = "msHidden";
    visibilityChange = "msvisibilitychange";
    } else if (typeof document.webkitHidden !== "undefined") {
    hidden = "webkitHidden";
    visibilityChange = "webkitvisibilitychange";
    }
    return {hidden, visibilityChange}
}

export function handleVisibilityChange(hidden) {
    let d = new Date().getTime();
    if (document[hidden]) {
        jatos.appendResultData(`Tab switch, component ${jatos.componentProperties.title}, ${d}\n`)
    } else {
        jatos.appendResultData(`Tab switch return, component ${jatos.componentProperties.title}, ${d}\n`)
    }
  }

export function suppressGoBack() {
    history.pushState(null, null, document.URL);
    window.addEventListener('popstate', function () {
        history.pushState(null, null, document.URL);
    });
}

export function getVect(n) {
    let vect = [];
    for (let i = 0; i < n; i++) {
        vect.push(i);
     }
     return vect;
}

export function shuffle(array) {
    for (let i = array.length - 1; i > 0; i--) {
      let j = Math.floor(Math.random() * (i + 1));
      [array[i], array[j]] = [array[j], array[i]];
    }
    return array;
}

export function indexOfAll(array, searchItem) {
    let i = array.indexOf(searchItem),
        indexes = [];
    while (i !== -1) {
      indexes.push(i);
      i = array.indexOf(searchItem, ++i);
    }
    return indexes;
  }

export function getRandomSubarray(arr, size) {
    let shuffled = arr.slice(0), i = arr.length, min = i - size, temp, index;
    while (i-- > min) {
        index = Math.floor((i + 1) * Math.random());
        temp = shuffled[index];
        shuffled[index] = shuffled[i];
        shuffled[i] = temp;
    }
    return shuffled.slice(min);
}  

export function JSON2CSV(objArray) {
  var array = typeof objArray != 'object' ? JSON.parse(objArray) : objArray;
  var str = '';

  for (var i = 0; i < array.length; i++) {
      var line = '';
      for (var index in array[i]) {
          if (line != '') line += ','

          line += array[i][index];
      }

      str += line + '\r\n';
  }

  return str;
}

export function getPics2NodesObj(vect) {
    let pics2nodesObj = {
        1: vect.slice(0, 11),
        // 2: vect.slice(14)
    }
    pics2nodesObj[1].unshift(null);
    // pics2nodesObj[2].unshift(null);
    return pics2nodesObj;
  }

export function isValid(value) {
    return typeof value;
}