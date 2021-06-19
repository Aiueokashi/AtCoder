function main(input){
    const [N,Q] = input[0].split(" ");
    const A = input[1].split(" ").map(i => BigInt(i))
    let C = new Array(N).fill(0);
    let p;
    for(let i = 0; i < N; i++){
         C[i] = A[i]-BigInt(i+1)
    }
    for(let i = 0; i < Q; i++){
        const K = BigInt(input[i+2]);
        const ans = iterativeFunction(C,K);
        if(ans === parseInt(N)){
            console.int(A[N-1] + (K - C[N-1]))
        }else{
            console.int((A[ans]-1n)-(C[ans]-K))
        }
    }
}
function iterativeFunction(arr, x) {
  let no = -1;
  let re = arr.length;
 
  while (re - no > 1) {
    const mid = no + Math.floor((re - no) / 2);
    if (arr[mid] >= x) re = mid;
    else no = mid;
  }
 
  return re;
}

console.int = function(bigint){
    console.log(bigint.toString())
}

main(require('fs').readFileSync('/dev/stdin', 'utf8').split("\n"));
