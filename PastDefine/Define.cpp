#define QLUBinarySearch(type,vec,k)  sort(all(vec)); \
int key = k; \
auto iter = type(all(vec),key);\
if (iter != vec.end() && *iter == key) {\
        cout << key << " is found at " << distance(vec.begin(), iter) << "\n";\
    } else {\
        cout << key << "is NOT found." << "\n";\
    }
    
#define QBinarySearch(vec, fkey, ekey) sort(all(vec));\
for (int key = fkey; key <= ekey; ++key) {\
        if (binary_search(vec.begin(), vec.end(), key)) {\
           print(key << ":" << "true") \
        } else {\
            print(key << ":" << "false")\
        }\
    }
