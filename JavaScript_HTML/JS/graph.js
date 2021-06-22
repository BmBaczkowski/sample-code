
export async function getAdjacencyMat() {
    let A = false;
    let str = `./matlabscripts/Adjacency.json`;
    let response = await fetch(str);
    if (response.ok === false) {
        return false;
    }
    A = await response.json();
    return A;
}

export async function getDistanceMat() {
    let A = false;
    let str = `./matlabscripts/Distance.json`;
    let response = await fetch(str);
    if (response.ok === false) {
        return false;
    }
    A = await response.json();
    return A;
}
