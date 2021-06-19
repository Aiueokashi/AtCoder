main(input){
    
}


try {
    main(require('fs').readFileSync('/dev/stdin', 'utf8'));
} catch (e) {
    try {
        main(require('fs').readFileSync('./dev/stdin', 'utf8'));
    } catch (e2) {
        
    }
}
