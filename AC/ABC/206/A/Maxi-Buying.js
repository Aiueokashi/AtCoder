function main(input){
    const N = Number(input[0])
    const c = N*1.08
    console.log(Math.trunc(c) < 206 ? "Yay!" : Math.trunc(c) === 206 ? "so-so" : ":(")
}

main(require('fs').readFileSync('/dev/stdin', 'utf8').split("\n"));
