#ifndef DATASTRUCTURE_HPP
#define DATASTRUCTURE_HPP
#include <vector>
#include <list>
#include <unordered_map>
#include <iterator>
#include <limits.h>
#include <string>
#include <tuple>

using namespace std;

struct HardBlock {
    std::tuple<int,int> block;
    int id;
    string name;
    HardBlock()
        : block(std::make_tuple(0, 0)),  // or std::tuple<int,int>(0,0)
          id(-1)
    {}
};


class Pad {
public:
    Pad() : id(-1), x(0), y(0) {};
    Pad(int id, int x, int y) : id(id), x(x), y(y) {};
    int id;
    int x;
    int y;
};

class Net {
public:
    Net() : id(-1), pad(nullptr) {};
    int id;
    Pad* pad;
    vector<HardBlock*> hardblocks;
};




#endif