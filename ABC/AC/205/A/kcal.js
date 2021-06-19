const input = require('fs').readFileSync('/dev/stdin', 'utf8');
const [n,m] = input.split(' ');
const _1cal = Number(n/100);
const ans = _1cal * Number(m);
console.log(ans)
