function main(input) {
    const [L,R] = input[0].split(" ").map(i => BigInt(i));
    let me = new Array(1000001).fill(0n);
    for (let i = 2; BigInt(i) <= R; ++i) {
        if (me[i] !== 0n)
            continue;
        for (let j = i; BigInt(j) <= R; j += i)
            me[j]++;
        if (BigInt(i) * BigInt(i) > BigInt(R))
            continue;
        for (let j = BigInt(i) * BigInt(i); j <= R; j += BigInt(i) * BigInt(i))
            me[j] = -R;
    }
    let count = 0n;
    for (let g = 2n; 2n * g <= R; ++g) {
        let cg = me[g];
        if (cg <= 0n)
            continue;
        let n1 = R / g - (L - 1n) / g;
        let n2 = n1 * (n1 - 1n) / 2n;
        count += (cg % 2n == 0n) ? -n2 : n2;
    }
    for (let g = 2n < L ? L : 2n ; 2n * g <= R; ++g) {
        count -= R / g - 1n;
    }
    count *= 2n
    console.log(count.toString())
}


main(require('fs').readFileSync('/dev/stdin', 'utf8').split("\n"));
