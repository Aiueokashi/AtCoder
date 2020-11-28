/*int main() {
  int N, A;
  cin >> N >> A;

  rep(i, N){
  string op;
  int B;
    int NUM;
    NUM = i+1;
    cin >> op >> B;
    if(op == "+"){
      A += B;
    }
    else if(op == "-"){
      A -= B;
    }
    else if(op == "*"){
      A *= B;
    }
    else if(op == "/"&&B != 0){
      A /= B;
    }
    else{
      cout << "error" << endl;
      break;
    }

    cout << NUM << ":" << A << endl;
    
  }
}*/