#include <fstream>
#include <iostream>
#include <vector>
#include <tuple>
#include <map>
#include <sstream>
#include <string>
#include <cstdlib>   // atoi
#include <algorithm>

using std::vector;
using std::tuple;
using std::get;
using std::string;
using std::map;
using std::ifstream;
using std::ofstream;
using std::cout;
using std::endl;

// Simple stoi replacement (g++ 4.8 libstdc++ is picky)
static int my_stoi(const string &s) {
    return std::atoi(s.c_str());
}

struct HardBlock {
    tuple<int,int> block; // (w, h)
    int id;
    string name;
    HardBlock() : block(std::make_tuple(0,0)), id(-1), name("") {}
};


// --------------------------------------------------
// ReadBlocks: format like ../testcase/sample2.txt:
//
// num_of_hardblocks 8
//
// // name w h
// hb0 4 2
// hb1 1 3
// ...
// --------------------------------------------------
bool ReadBlocks(const string &path, vector<HardBlock> &blocks) {
    ifstream fin(path.c_str());
    if (!fin) {
        std::cerr << "Cannot open blocks file: " << path << std::endl;
        return false;
    }

    blocks.clear();

    string token;
    int N = 0;

    // Find "num_of_hardblocks N"
    while (fin >> token) {
        if (token[0] == '/' || token[0] == '#') {
            // comment line, skip rest of line
            std::string dummy;
            std::getline(fin, dummy);
            continue;
        }
        if (token == "num_of_hardblocks") {
            fin >> N;
            break;
        } else {
            // Unexpected token; ignore this line
            std::string dummy;
            std::getline(fin, dummy);
        }
    }

    if (N <= 0) {
        std::cerr << "Failed to read num_of_hardblocks from " << path << std::endl;
        return false;
    }

    // Now read N lines with: hbX w h (skip comments/blank lines)
    int count = 0;
    while (count < N && fin.good()) {
        std::string line;
        if (!std::getline(fin, line)) break;
        if (line.empty()) continue;
        if (line[0] == '/' || line[0] == '#') continue;

        std::istringstream iss(line);
        string name;
        int w, h;
        if (!(iss >> name >> w >> h)) {
            // malformed line, skip
            continue;
        }
        HardBlock hb;
        hb.block = std::make_tuple(w, h);
        hb.id    = count;
        hb.name  = name;   // <-- store "hb0", "hb1", ...
        blocks.push_back(hb);
        ++count;
    }

    if ((int)blocks.size() != N) {
        std::cerr << "Warning: expected " << N << " blocks, but read "
                  << blocks.size() << " from " << path << std::endl;
    }

    std::cerr << "ReadBlocks: " << blocks.size() << " blocks\n";
    return !blocks.empty();
}

// --------------------------------------------------
// ReadSolution: your sample2.sol format:
//
// # Objective value = 7
// y 7
// x1 6
// r1 0
// y1 0
// ...
// --------------------------------------------------
bool ReadSolution(const string &path,
                  double &chipHeight,
                  vector<double> &x,
                  vector<double> &y,
                  vector<int> &r)
{
    ifstream fin(path.c_str());
    if (!fin) {
        std::cerr << "Cannot open solution file: " << path << std::endl;
        return false;
    }

    chipHeight = 0.0;
    map<int,double> xMap;
    map<int,double> yMap;
    map<int,int>    rMap;

    int N = 0;
    string line;

    while (std::getline(fin, line)) {
        if (line.empty()) continue;
        if (line[0] == '#') continue;

        std::istringstream iss(line);
        string var;
        double val;
        if (!(iss >> var >> val)) continue;

        if (var == "y") {
            chipHeight = val;
            continue;
        }

        char c = var[0];
        if (c == 'x' || c == 'y' || c == 'r') {
            string idxStr = var.substr(1);
            int idx = my_stoi(idxStr);  // 1-based
            if (idx <= 0) continue;
            if (idx > N) N = idx;

            if (c == 'x')      xMap[idx] = val;
            else if (c == 'y') yMap[idx] = val;
            else if (c == 'r') rMap[idx] = (val != 0.0) ? 1 : 0;
        }
    }

    if (N == 0) {
        std::cerr << "No x_i/y_i/r_i found in solution file\n";
        return false;
    }

    x.assign(N, 0.0);
    y.assign(N, 0.0);
    r.assign(N, 0);

    for (map<int,double>::const_iterator it = xMap.begin(); it != xMap.end(); ++it) {
        int idx = it->first;
        if (idx >= 1 && idx <= N) x[idx-1] = it->second;
    }
    for (map<int,double>::const_iterator it = yMap.begin(); it != yMap.end(); ++it) {
        int idx = it->first;
        if (idx >= 1 && idx <= N) y[idx-1] = it->second;
    }
    for (map<int,int>::const_iterator it = rMap.begin(); it != rMap.end(); ++it) {
        int idx = it->first;
        if (idx >= 1 && idx <= N) r[idx-1] = it->second;
    }

    std::cerr << "ReadSolution: N=" << N << ", chipHeight=" << chipHeight << "\n";
    return true;
}

// --------------------------------------------------
// VisualizeFloorplanSVG
// --------------------------------------------------
void VisualizeFloorplanSVG(const vector<HardBlock>& blocks,
                           const vector<double>& x,
                           const vector<double>& y,
                           const vector<int>& r,
                           double chipHeight,
                           const string& svg_path)
{
    size_t N = blocks.size();
    if (x.size() != N || y.size() != N || r.size() != N) {
        std::cerr << "Size mismatch: blocks=" << N
                  << " x=" << x.size()
                  << " y=" << y.size()
                  << " r=" << r.size() << std::endl;
        return;
    }
    if (chipHeight <= 0) {
        std::cerr << "Invalid chipHeight=" << chipHeight << std::endl;
        return;
    }

    // compute chip width = max(x_i + width_i)
    double chipWidth = 0.0;
    for (size_t i = 0; i < N; ++i) {
        int wi = get<0>(blocks[i].block);
        int hi = get<1>(blocks[i].block);
        double bw = (r[i] == 0) ? wi : hi;
        double candidate = x[i] + bw;
        if (candidate > chipWidth) chipWidth = candidate;
    }
    if (chipWidth <= 0) chipWidth = 1.0;

    double maxDim = (chipWidth > chipHeight) ? chipWidth : chipHeight;
    double targetPixels = 1000.0;
    double scale = targetPixels / maxDim;
    if (scale <= 0) scale = 1.0;

    double svgW = chipWidth  * scale;
    double svgH = chipHeight * scale;

    ofstream ofs(svg_path.c_str());
    if (!ofs) {
        std::cerr << "Cannot open SVG for write: " << svg_path << std::endl;
        return;
    }

    ofs << "<?xml version=\"1.0\" standalone=\"no\"?>\n";
    ofs << "<svg xmlns=\"http://www.w3.org/2000/svg\" version=\"1.1\""
        << " width=\"" << svgW << "\""
        << " height=\"" << svgH << "\""
        << " viewBox=\"0 0 " << svgW << " " << svgH << "\">\n\n";

    ofs << "  <rect x=\"0\" y=\"0\" width=\"" << svgW
        << "\" height=\"" << svgH
        << "\" fill=\"none\" stroke=\"black\" stroke-width=\"1\" />\n\n";

    const char* colors[6] = {
        "#ffcccc", "#ccffcc", "#ccccff",
        "#ffe0b3", "#e0ccff", "#ccffff"
    };

    for (size_t i = 0; i < N; ++i) {
        int wi = get<0>(blocks[i].block);
        int hi = get<1>(blocks[i].block);

        double bw = (r[i] == 0) ? wi : hi;
        double bh = (r[i] == 0) ? hi : wi;

        double x_model = x[i];
        double y_model = y[i];

        double x_svg = x_model * scale;
        double y_svg = (chipHeight - (y_model + bh)) * scale;
        double w_svg = bw * scale;
        double h_svg = bh * scale;

        const char* fillColor = colors[i % 6];

        ofs << "  <rect x=\"" << x_svg
            << "\" y=\"" << y_svg
            << "\" width=\"" << w_svg
            << "\" height=\"" << h_svg
            << "\" fill=\"" << fillColor
            << "\" stroke=\"black\" stroke-width=\"0.5\" />\n";

        // Label blocks with their ID/name (e.g., "hb0") for moderate N
        if (N <= 200) {
            double cx = x_svg + w_svg * 0.5;
            double cy = y_svg + h_svg * 0.5;

            // Prefer the name (hb0, hb1, ...) if available
            string label = blocks[i].name;
            if (label.empty()) {
                std::ostringstream oss;
                oss << (i + 1);
                label = oss.str();
            }

            // font size in SVG units (roughly pixels)
            double fontSize = std::min(svgW, svgH) / 30.0; // ~33px if 1000x1000
            if (fontSize < 10.0) fontSize = 10.0;          // minimum readable

            ofs << "  <text x=\"" << cx
                << "\" y=\"" << cy
                << "\" font-family=\"sans-serif\""
                << " font-size=\"" << fontSize << "\""
                << " fill=\"black\""
                << " text-anchor=\"middle\""
                << " dominant-baseline=\"middle\">"
                << label << "</text>\n";
        }
    }

    ofs << "\n</svg>\n";
    ofs.close();

    std::cerr << "Floorplan SVG written to " << svg_path << std::endl;
}

// --------------------------------------------------
// main
// --------------------------------------------------
int main(int argc, char** argv) {
    if (argc < 4) {
        std::cerr << "Usage: " << argv[0]
                  << " <blocks_file> <solution_file> <output_svg>\n";
        return 1;
    }

    string blockPath = argv[1];
    string solPath   = argv[2];
    string svgPath   = argv[3];

    vector<HardBlock> blocks;
    if (!ReadBlocks(blockPath, blocks)) {
        return 1;
    }

    double chipHeight = 0.0;
    vector<double> x, y;
    vector<int> r;
    if (!ReadSolution(solPath, chipHeight, x, y, r)) {
        return 1;
    }

    if (blocks.size() != x.size()) {
        std::cerr << "Warning: blocks.size()=" << blocks.size()
                  << " != x.size()=" << x.size() << std::endl;
    }

    VisualizeFloorplanSVG(blocks, x, y, r, chipHeight, svgPath);
    return 0;
}
