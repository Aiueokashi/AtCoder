function main(input){
    const N = parseInt(input[0]);
    if(N <= 2){
        return console.log(N)
    }
    for(let i = 1; i < N; i++){
        if(i*(i+1)/2 >= N){
            return console.log(i)
        }
    }
}

main(require('fs').readFileSync('/dev/stdin', 'utf8').split("\n"));
