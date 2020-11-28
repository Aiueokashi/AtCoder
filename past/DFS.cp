ll N;
int dfs(string &S) {
    ll ret = 0;
    vector<char> lis = {'3', '5', '7'};
    if ((int)S.size() > 0) {
        if (stoll(S) > N) {   
            return ret;
        } else {  
            bool ok = true;
            for (auto c : lis) {
                if (S.find(c) == string::npos) ok = false;
            }
            if (ok) ret++;
        }
    }
    for (auto c : lis) {
        S.push_back(c);
        ret += dfs(S);
        S.pop_back();
    }
    return ret;
}
