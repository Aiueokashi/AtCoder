const input = require('fs').readFileSync('/dev/stdin', 'utf8').split('\n');
const num = parseInt(input[0]);
const [...array] = input[1].split(' ');
let result = [];
array.forEach(a => {
	result[Number(a)] = true;
})
//console.log(result)
for(let i = 0; i < num; i++){
    if(!result[i+1]){
        console.log("No")
        return
    }
}
console.log("Yes");
