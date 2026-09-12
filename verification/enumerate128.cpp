#include <cstdint>
#include <fstream>
#include <iostream>
#include <vector>
using U=unsigned __int128;
int main(int argc,char**argv){std::ifstream f(argv[1]);int n,k;f>>n>>k;if(n>128||k>28||k<1)return 2;std::vector<U>b(k);for(auto&v:b){uint64_t lo,hi;f>>lo>>hi;v=U(lo)|(U(hi)<<64);}U v=0,bestv=0;int best=n;uint64_t count=0,total=1ULL<<k;for(uint64_t i=1;i<total;i++){v^=b[__builtin_ctzll(i)];int w=__builtin_popcountll((uint64_t)v)+__builtin_popcountll((uint64_t)(v>>64));if(w<best){best=w;bestv=v;count=1;}else if(w==best)count++;}std::cout<<"{\"complete\":true,\"k\":"<<k<<",\"n\":"<<n<<",\"d\":"<<best<<",\"count\":"<<count<<",\"nonzero_words_checked\":"<<total-1<<",\"support\":[";bool sep=false;for(int t=0;t<n;t++)if((bestv>>t)&1){if(sep)std::cout<<",";std::cout<<t;sep=true;}std::cout<<"]}\n";}
