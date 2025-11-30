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

struct Net {
    string name;               // e.g. "net0", "nDIR"
    vector<int> hb_indices;    // indices into HardBlock array (0-based)
    vector<int> pad_indices;   // pad indices (0..num_of_pads-1)
};

static string trim(const string &s) {
    size_t start = s.find_first_not_of(" \t\r\n");
    if (start == string::npos) return "";
    size_t end = s.find_last_not_of(" \t\r\n");
    return s.substr(start, end - start + 1);
}

// --------------------------------------------------
// ReadBlocks: format like ami33_WL.txt:
//
// num_of_hardblocks 33
// num_of_nets 43
// num_of_pads 4
//
// // hardblocks
// // name w h
// hbk1   413 220
// hbk2   112 321
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

    // Now read lines with: <name> w h  (hbk1 413 220, hb0 4 2, etc.)
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
        hb.name  = name;   // e.g. "hbk1" or "hb0"
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
// ReadNets: from the same testcase file, after headers.
//
// num_of_hardblocks 33
// num_of_nets 43
// num_of_pads 4
// ...
// hbk1 413 220
// ...
// nDIR hbk9a hbk11
// nC30 hbk8a hbk8b
// ...
// --------------------------------------------------
bool ReadNets(const string &path,
              const vector<HardBlock> &blocks,
              vector<Net> &nets,
              int &num_of_pads)
{
    ifstream fin(path.c_str());
    if (!fin) {
        std::cerr << "Cannot open nets file: " << path << std::endl;
        return false;
    }

    nets.clear();
    num_of_pads = 0;

    // For mapping "<blockName>" -> index
    map<string,int> hbNameToIdx;
    for (size_t i = 0; i < blocks.size(); ++i) {
        hbNameToIdx[blocks[i].name] = (int)i;
    }

    string line;
    int num_of_nets_expected = 0;

    while (std::getline(fin, line)) {
        string t = trim(line);
        if (t.empty()) continue;
        if (t[0] == '/' || t[0] == '#') continue;

        std::istringstream iss(t);
        string first;
        iss >> first;
        if (first.empty()) continue;

        if (first == "num_of_nets") {
            iss >> num_of_nets_expected;
            continue;
        }
        if (first == "num_of_pads") {
            iss >> num_of_pads;
            continue;
        }
        if (first == "num_of_hardblocks") {
            // header line, already handled in ReadBlocks
            continue;
        }

        // If the first token is a block name, it's a block definition line, not a net.
        if (hbNameToIdx.find(first) != hbNameToIdx.end()) {
            continue;
        }

        // Otherwise treat this line as a net definition:
        // e.g. "net0 hb0 hb1 pad0", "nDIR hbk9a hbk11", etc.
        Net net;
        net.name = first;  // store raw net name

        string pin;
        while (iss >> pin) {
            // Is it a block pin? (match by name)
            auto it = hbNameToIdx.find(pin);
            if (it != hbNameToIdx.end()) {
                net.hb_indices.push_back(it->second);
                continue;
            }

            // Is it a pad pin? "pad0".."pad3"
            if (pin.rfind("pad", 0) == 0) {
                string pStr = pin.substr(3);
                int pidx = my_stoi(pStr);
                net.pad_indices.push_back(pidx);
                continue;
            }

            // Common power/ground nets we want to ignore
            if (pin == "GND" || pin == "PWR" || pin == "VDD" || pin == "VSS") {
                continue;
            }

            // Otherwise, unknown token in this net line
            std::cerr << "Warning: unknown net pin '" << pin
                      << "' in line: " << t << std::endl;
        }

        nets.push_back(net);
    }

    if (num_of_nets_expected > 0 && (int)nets.size() != num_of_nets_expected) {
        std::cerr << "Warning: expected " << num_of_nets_expected
                  << " nets, but read " << nets.size() << " from " << path << std::endl;
    }

    std::cerr << "ReadNets: " << nets.size()
              << " nets, num_of_pads=" << num_of_pads << "\n";
    return true;
}

// --------------------------------------------------
// ReadSolution: .sol format:
//
// # Objective value = 7
// y 7
// x1 6
// r1 0
// y1 0
// ...
// plus many extra vars (max_...x, min_...y, p12, q12, ...)
// we only care about x_i, y_i, r_i and global y.
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
        if (line[0] == '#') continue;  // comment

        std::istringstream iss(line);
        string var;
        double val;
        if (!(iss >> var >> val)) continue;

        // global chip height
        if (var == "y") {
            chipHeight = val;
            continue;
        }

        // We only care about x<idx>, y<idx>, r<idx>
        if (var.size() < 2) continue;

        char c = var[0];
        if (c != 'x' && c != 'y' && c != 'r') {
            // e.g. max_nC30x, min_nC30y, p12, q12 -> ignore
            continue;
        }

        string idxStr = var.substr(1);
        if (idxStr.empty() || !std::isdigit(static_cast<unsigned char>(idxStr[0]))) {
            continue;
        }

        int idx = my_stoi(idxStr);  // 1-based
        if (idx <= 0) continue;
        if (idx > N) N = idx;

        if (c == 'x')      xMap[idx] = val;
        else if (c == 'y') yMap[idx] = val;
        else if (c == 'r') rMap[idx] = (val != 0.0) ? 1 : 0;
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

void VisualizeFloorplanSVG(const vector<HardBlock>& blocks,
                           const vector<double>& x,
                           const vector<double>& y,
                           const vector<int>& r,
                           double chipHeight,
                           const vector<Net>& nets,
                           int num_of_pads,
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

    // chip border
    ofs << "  <rect x=\"0\" y=\"0\" width=\"" << svgW
        << "\" height=\"" << svgH
        << "\" fill=\"none\" stroke=\"black\" stroke-width=\"1\" />\n\n";

    const char* blockColors[6] = {
        "#ffcccc", "#ccffcc", "#ccccff",
        "#ffe0b3", "#e0ccff", "#ccffff"
    };

    // draw blocks
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

        const char* fillColor = blockColors[i % 6];

        ofs << "  <rect x=\"" << x_svg
            << "\" y=\"" << y_svg
            << "\" width=\"" << w_svg
            << "\" height=\"" << h_svg
            << "\" fill=\"" << fillColor
            << "\" stroke=\"black\" stroke-width=\"0.5\" />\n";

        // Label blocks with their name (hbk1, hb0, ...) for moderate N
        if (N <= 200) {
            double cx = x_svg + w_svg * 0.5;
            double cy = y_svg + h_svg * 0.5;

            string label = blocks[i].name;
            if (label.empty()) {
                std::ostringstream oss;
                oss << (i + 1);
                label = oss.str();
            }

            double fontSize = std::min(svgW, svgH) / 30.0;
            if (fontSize < 10.0) fontSize = 10.0;

            if (r[i] == 0) {
                // normal (not rotated): horizontal text
                ofs << "  <text x=\"" << cx
                    << "\" y=\"" << cy
                    << "\" font-family=\"sans-serif\""
                    << " font-size=\"" << fontSize << "\""
                    << " fill=\"black\""
                    << " text-anchor=\"middle\""
                    << " dominant-baseline=\"middle\">"
                    << label << "</text>\n";
            } else {
                // rotated block: draw label rotated -90 degrees around its center
                ofs << "  <text x=\"" << cx
                    << "\" y=\"" << cy
                    << "\" font-family=\"sans-serif\""
                    << " font-size=\"" << fontSize << "\""
                    << " fill=\"black\""
                    << " text-anchor=\"middle\""
                    << " dominant-baseline=\"middle\""
                    << " transform=\"rotate(-90 " << cx << " " << cy << ")\">"
                    << label << "</text>\n";
            }
        }
    }


    ofs << "\n";

    // --------------------------------------------------------
    // Pre-compute pad positions in *model* coords (middle of edges)
    // To match LP generator mapping:
    // pad0: left edge middle   -> (0, floor(chipHeight/2))
    // pad1: top edge middle    -> (floor(chipWidth/2), chipHeight)
    // pad2: right edge middle  -> (chipWidth, floor(chipHeight/2))
    // pad3: bottom edge middle -> (floor(chipWidth/2), 0)
    // --------------------------------------------------------
    double midX_int = (double)((int)(chipWidth  / 2.0));
    double midY_int = (double)((int)(chipHeight / 2.0));

    struct PadPos { double x, y; };
    vector<PadPos> padPos(num_of_pads);

    if (num_of_pads >= 1) padPos[0] = { 0.0,      midY_int };
    if (num_of_pads >= 2) padPos[1] = { midX_int, chipHeight };
    if (num_of_pads >= 3) padPos[2] = { chipWidth, midY_int };
    if (num_of_pads >= 4) padPos[3] = { midX_int, 0.0 };

    // --------------------------------------------------------
    // Draw nets: connect each pin to gravity center (no circles)
    // Each net uses a different color.
    // --------------------------------------------------------
    const char* netColors[] = {
        "#0000cc", "#cc0000", "#009900",
        "#9933cc", "#ff9900", "#0099cc",
        "#cc0099", "#666600"
    };
    const int numNetColors = sizeof(netColors) / sizeof(netColors[0]);

    for (size_t ni = 0; ni < nets.size(); ++ni) {
        const Net &net = nets[ni];

        // collect all pin coordinates (model space)
        vector<double> pxs;
        vector<double> pys;

        // block pins: center of placed rectangle
        for (size_t k = 0; k < net.hb_indices.size(); ++k) {
            int idx = net.hb_indices[k];
            if (idx < 0 || idx >= (int)N) continue;

            int wi = get<0>(blocks[idx].block);
            int hi = get<1>(blocks[idx].block);
            double bw = (r[idx] == 0) ? wi : hi;
            double bh = (r[idx] == 0) ? hi : wi;

            double cx = x[idx] + bw * 0.5;
            double cy = y[idx] + bh * 0.5;

            pxs.push_back(cx);
            pys.push_back(cy);
        }

        // pad pins: precomputed padPos
        for (size_t k = 0; k < net.pad_indices.size(); ++k) {
            int pid = net.pad_indices[k];
            if (pid < 0 || pid >= num_of_pads) continue;
            pxs.push_back(padPos[pid].x);
            pys.push_back(padPos[pid].y);
        }

        if (pxs.empty()) continue;

        // gravity center (centroid) of all pins
        double sumx = 0.0, sumy = 0.0;
        for (size_t i = 0; i < pxs.size(); ++i) {
            sumx += pxs[i];
            sumy += pys[i];
        }
        double cx = sumx / pxs.size();
        double cy = sumy / pys.size();

        double cx_svg = cx * scale;
        double cy_svg = (chipHeight - cy) * scale;

        // choose net color
        string netColor = netColors[ni % numNetColors];

        // draw lines from centroid to each pin (slightly thicker)
        for (size_t i = 0; i < pxs.size(); ++i) {
            double px_svg = pxs[i] * scale;
            double py_svg = (chipHeight - pys[i]) * scale;

            ofs << "  <line x1=\"" << cx_svg
                << "\" y1=\"" << cy_svg
                << "\" x2=\"" << px_svg
                << "\" y2=\"" << py_svg
                << "\" stroke=\"" << netColor
                << "\" stroke-width=\"2\" stroke-opacity=\"0.7\" />\n";
        }
    }

    // --------------------------------------------------------
    // Draw pad markers and labels on top so they're visible
    // (label directions adjusted to match new mapping)
    // --------------------------------------------------------
    for (int pid = 0; pid < num_of_pads; ++pid) {
        double px = padPos[pid].x;
        double py = padPos[pid].y;
        double px_svg = px * scale;
        double py_svg = (chipHeight - py) * scale;

        double padSize = 6.0; // in SVG units (~pixels)
        double fontSize = 12.0;

        ofs << "  <rect x=\"" << (px_svg - padSize/2.0)
            << "\" y=\"" << (py_svg - padSize/2.0)
            << "\" width=\"" << padSize
            << "\" height=\"" << padSize
            << "\" fill=\"#ff0000\" stroke=\"black\" stroke-width=\"0.5\" />\n";

        double textX = px_svg;
        double textY = py_svg;

        // Place label depending on which edge the pad is on
        if (pid == 0) {
            // left middle: label to the right
            textX = px_svg + padSize + 2;
            textY = py_svg + fontSize / 2.0;
        } else if (pid == 1) {
            // top center: label below
            textY = py_svg + padSize + fontSize;
        } else if (pid == 2) {
            // right middle: label to the left
            textX = px_svg - padSize - 25;
            textY = py_svg + fontSize / 2.0;
        } else if (pid == 3) {
            // bottom center: label above
            textY = py_svg - padSize - 2;
        }

        ofs << "  <text x=\"" << textX
            << "\" y=\"" << textY
            << "\" font-family=\"sans-serif\""
            << " font-size=\"" << fontSize << "\""
            << " fill=\"black\">"
            << "pad" << pid << "</text>\n";
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

    // Read nets + pads from same testcase file
    vector<Net> nets;
    int num_of_pads = 0;
    if (!ReadNets(blockPath, blocks, nets, num_of_pads)) {
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

    VisualizeFloorplanSVG(blocks, x, y, r, chipHeight, nets, num_of_pads, svgPath);
    return 0;
}
