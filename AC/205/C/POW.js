const input = require('fs').readFileSync('/dev/stdin', 'utf8').split(' ');

for(let i = 0; i < input.length; i ++){
    eval(`var ${String.fromCharCode(i+65)} = ${parseInt(input[i])}`)
}
let arr = ["A","B","C"]
for(let i = 0; i < arr.length; i++){
    eval(`var ${arr[i].toLowerCase()} = ${Math.abs(parseInt(input[i]))}`)
}
if(C%2 === 0){
    console.log(a === b ? "=" : a > b ? ">" : "<")
}else{
    console.log(A === B ? "=" : A > B ? ">" : "<")
}
