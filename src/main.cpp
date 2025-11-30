#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <tuple>
#include <ctime>
#include <sys/time.h>
#include <time.h>
#include <cmath>
#include <random>
#include <algorithm>
#include <cctype>       // std::isdigit, std::isalpha, ...
#include <unordered_map>

#include "datastructure.hpp"

#define EQ 0
#define GE 1
#define LE 2

using namespace std;

// -------------------------
// Helper structs for HPWL
// -------------------------
struct NetHPWL {
    std::string name;        // 原始 net 名字，例如 net0, nDIR, nC30, ...
    std::vector<int> block_ids; // indices into hardblocks array
    std::vector<int> pad_ids;   // pad indices 0..num_of_pads-1
};

double GetTime(void) {
    struct timeval tv;
    gettimeofday(&tv, NULL);
    return tv.tv_sec + tv.tv_usec * 1e-6;
}

// map (i,j) with 0 <= i < j < N to [0, numPermu-1]
static int pairIndex(int i, int j, int N) {
    // sum_{k=0}^{i-1} (N-k-1) + (j-i-1)
    return i * (2 * N - i - 1) / 2 + (j - i - 1);
}

// extract numeric part from a string like "net5" -> 5
static int extractNetId(const std::string& name) {
    std::string digits;
    for (char c : name) {
        if (std::isdigit(static_cast<unsigned char>(c))) {
            digits.push_back(c);
        }
    }
    if (!digits.empty()) return std::stoi(digits);
    return -1; // fallback
}

// Make a safe tag for LP variable names from a net name.
// Special case: keep old behavior for "net<digits>" → "n<digits>".
static std::string makeNetTag(const std::string& name) {
    // keep backward compatible for net0/net1/... → n0/n1/...
    if (name.rfind("net", 0) == 0) {
        std::string digits;
        for (char c : name) {
            if (std::isdigit(static_cast<unsigned char>(c))) {
                digits.push_back(c);
            }
        }
        if (!digits.empty()) {
            return "n" + digits;  // old style: max_n0x, max_n1x, ...
        }
    }

    // general case: sanitize the name itself
    std::string tag;
    tag.reserve(name.size() + 1);

    // ensure first char is a letter or '_'
    if (name.empty() ||
        !(std::isalpha(static_cast<unsigned char>(name[0])) || name[0] == '_')) {
        tag.push_back('n');
    }

    for (char c : name) {
        if (std::isalnum(static_cast<unsigned char>(c)) || c == '_') {
            tag.push_back(c);
        } else {
            tag.push_back('_');
        }
    }
    return tag;
}

// -------------------------------------------------------
// Sutanthavibul: dimension + non-overlap -> matrix
// -------------------------------------------------------
void Sutanthavibul(HardBlock* hardblocks, int N, int W, vector<vector<int>>& matrix) {
    int numPermu = N * (N - 1) / 2;

    // #variables (without OP/RHS) = 3N + 2*numPermu + 1(global y)
    // = N^2 + 2N + 1 = (N+1)^2
    int number_of_variables = N * N + 2 * N + 3;  // +1 for global y, +OP +RHS
    int number_of_constraints = 2 * (N * N);      // 2N^2

    matrix.assign(number_of_constraints,
                  vector<int>(number_of_variables, 0));

    // column indices
    auto idx_x = [N](int i) { return i; };                     // x_i
    auto idx_y = [N](int i) { return N + i; };                 // y_i
    auto idx_r = [N](int i) { return 2 * N + i; };             // r_i

    auto idx_p = [N](int i, int j) {
        return 3 * N + pairIndex(i, j, N);                     // p_ij
    };
    auto idx_q = [N, numPermu](int i, int j) {
        return 3 * N + numPermu + pairIndex(i, j, N);          // q_ij
    };

    int y_global = 3 * N + 2 * numPermu;                       // y
    int OP       = number_of_variables - 2;                    // operator
    int RHS      = number_of_variables - 1;                    // RHS

    // Big-M
    int M = 0;
    for (int i = 0; i < N; ++i) {
        int wi = get<0>(hardblocks[i].block);
        int hi = get<1>(hardblocks[i].block);
        M += max(wi, hi);
    }
    std::cout << "M: " << M << std::endl;

    int row = 0;

    // ---------------------------------------------------
    // 1) Dimension constraints (x and y) for each block i
    // ---------------------------------------------------
    for (int i = 0; i < N; ++i) {
        int wi = get<0>(hardblocks[i].block);
        int hi = get<1>(hardblocks[i].block);

        // (1) x_i + r_i h_i + (1-r_i) w_i <= W
        //     -> x_i + r_i(h_i - w_i) <= W - w_i
        {
            int c_xi = idx_x(i);
            int c_ri = idx_r(i);

            matrix[row][c_xi] = 1;
            matrix[row][c_ri] = hi - wi;
            matrix[row][OP]   = LE;
            matrix[row][RHS]  = W - wi;
            ++row;
        }

        // (2) y_i + r_i w_i + (1-r_i) h_i <= y
        //     -> y_i - y + r_i(w_i - h_i) <= -h_i
        {
            int c_yi = idx_y(i);
            int c_ri = idx_r(i);

            matrix[row][c_yi]     = 1;
            matrix[row][y_global] = -1;          // -y
            matrix[row][c_ri]     = wi - hi;
            matrix[row][OP]       = LE;
            matrix[row][RHS]      = -hi;
            ++row;
        }
    }

    // ---------------------------------------------------
    // 2) Non-overlap constraints for each pair (i,j), i<j
    // ---------------------------------------------------
    for (int i = 0; i < N; ++i) {
        int wi = get<0>(hardblocks[i].block);
        int hi = get<1>(hardblocks[i].block);

        for (int j = i + 1; j < N; ++j) {
            int wj = get<0>(hardblocks[j].block);
            int hj = get<1>(hardblocks[j].block);

            int c_xi = idx_x(i);
            int c_xj = idx_x(j);
            int c_yi = idx_y(i);
            int c_yj = idx_y(j);
            int c_ri = idx_r(i);
            int c_rj = idx_r(j);
            int c_p  = idx_p(i, j);
            int c_q  = idx_q(i, j);

            // (3) x_i + r_i h_i + (1-r_i) w_i <= x_j + M(p_ij + q_ij)
            //     -> x_i - x_j + r_i(h_i - w_i) - M p_ij - M q_ij <= -w_i
            matrix[row][c_xi] =  1;
            matrix[row][c_xj] = -1;
            matrix[row][c_ri] =  hi - wi;
            matrix[row][c_p]  = -M;
            matrix[row][c_q]  = -M;
            matrix[row][OP]   = LE;
            matrix[row][RHS]  = -wi;
            ++row;

            // (4) y_i + r_i w_i + (1-r_i) h_i <= y_j + M(1 + p_ij - q_ij)
            //     -> y_i - y_j + r_i(w_i - h_i) - M p_ij + M q_ij <= M - h_i
            matrix[row][c_yi] =  1;
            matrix[row][c_yj] = -1;
            matrix[row][c_ri] =  wi - hi;
            matrix[row][c_p]  = -M;
            matrix[row][c_q]  =  M;
            matrix[row][OP]   = LE;
            matrix[row][RHS]  = M - hi;
            ++row;

            // (5) x_i - r_j h_j - (1-r_j) w_j >= x_j - M(1 - p_ij + q_ij)
            //     -> x_i - x_j + r_j(w_j - h_j) - M p_ij + M q_ij >= -M + w_j
            matrix[row][c_xi] =  1;
            matrix[row][c_xj] = -1;
            matrix[row][c_rj] =  wj - hj;
            matrix[row][c_p]  = -M;
            matrix[row][c_q]  =  M;
            matrix[row][OP]   = GE;
            matrix[row][RHS]  = -M + wj;
            ++row;

            // (6) y_i - r_j w_j - (1-r_j) h_j >= y_j - M(2 - p_ij - q_ij)
            //     -> y_i - y_j + r_j(h_j - w_j) - M p_ij - M q_ij >= -2M + hj
            matrix[row][c_yi] =  1;
            matrix[row][c_yj] = -1;
            matrix[row][c_rj] =  hj - wj;
            matrix[row][c_p]  = -M;
            matrix[row][c_q]  = -M;
            matrix[row][OP]   = GE;
            matrix[row][RHS]  = -2 * M + hj;
            ++row;
        }
    }
}

// -------------------------------------------------------
// Convert matrix + HPWL (blocks + pads) into LP file
// -------------------------------------------------------
void ConvertMatrixToLPFILE(vector<vector<int>>& matrix,
                           const std::string& output_filepath,
                           HardBlock* hardblocks,
                           int num_of_hardblocks,
                           const std::vector<NetHPWL>& nets_hpwl,
                           int dieWidth) {
    ofstream output_file(output_filepath);
    if (!output_file) return;
    if (matrix.empty()) {
        output_file.close();
        return;
    }

    int numRows = (int)matrix.size();
    int numCols = (int)matrix[0].size();

    int RHScol = numCols - 1;
    int OPcol  = numCols - 2;
    int ycol   = numCols - 3; // global y

    // deduce N from #constraints: 2N^2 = numRows
    int N = (int)std::lround(std::sqrt(numRows / 2.0));
    int numPermu = N * (N - 1) / 2;

    if (N != num_of_hardblocks) {
        std::cerr << "Warning: deduced N = " << N
                  << " but num_of_hardblocks = " << num_of_hardblocks << "\n";
    }

    auto idx_x = [N](int i) { return i; };
    auto idx_y = [N](int i) { return N + i; };
    auto idx_r = [N](int i) { return 2 * N + i; };

    auto idx_p = [N](int i, int j) {
        return 3 * N + pairIndex(i, j, N);
    };
    auto idx_q = [N, numPermu](int i, int j) {
        return 3 * N + numPermu + pairIndex(i, j, N);
    };

    // build variable names for all columns
    vector<string> varName(numCols, "");

    for (int i = 0; i < N; ++i) {
        varName[idx_x(i)] = "x" + std::to_string(i + 1);
        varName[idx_y(i)] = "y" + std::to_string(i + 1);
        varName[idx_r(i)] = "r" + std::to_string(i + 1);
    }

    for (int i = 0; i < N; ++i) {
        for (int j = i + 1; j < N; ++j) {
            int cp = idx_p(i, j);
            int cq = idx_q(i, j);
            varName[cp] = "p" + std::to_string(i + 1) + std::to_string(j + 1);
            varName[cq] = "q" + std::to_string(i + 1) + std::to_string(j + 1);
        }
    }

    varName[ycol] = "y"; // global chip height

    // ---------- Minimize ----------
    output_file << "Minimize\n";
    output_file << "  20 y";

    // Add Σ_n (max_nx - min_nx + max_ny - min_ny)
    for (const auto& net : nets_hpwl) {
        std::string tag = makeNetTag(net.name);
        output_file << " + max_" << tag << "x"
                    << " - min_" << tag << "x"
                    << " + max_" << tag << "y"
                    << " - min_" << tag << "y";
    }
    output_file << "\n";

    // ---------- Subject To ----------
    output_file << "Subject To\n";

    // 1) Existing Sutanthavibul constraints from matrix
    for (int r = 0; r < numRows; ++r) {
        output_file << "  c" << (r + 1) << ": ";

        bool firstTerm = true;

        for (int c = 0; c < numCols; ++c) {
            if (c == OPcol || c == RHScol) continue;

            int coef = matrix[r][c];
            if (coef == 0) continue;
            const string& name = varName[c];
            if (name.empty()) continue;

            int abscoef = std::abs(coef);

            if (firstTerm) {
                if (coef < 0) output_file << "- ";
            } else {
                output_file << (coef < 0 ? " - " : " + ");
            }
            if (abscoef != 1) {
                output_file << abscoef << " ";
            }
            output_file << name;
            firstTerm = false;
        }

        if (firstTerm) {
            // no variable term in this row
            output_file << "0";
        }

        int op = matrix[r][OPcol];
        if (op == LE)      output_file << " <= ";
        else if (op == GE) output_file << " >= ";
        else               output_file << " = ";

        output_file << matrix[r][RHScol] << "\n";
    }

    // 2) HPWL constraints: for each net, for each pin (HB or pad)
    int cid = numRows + 1;

    double W      = static_cast<double>(dieWidth);
    double W_mid  = std::floor(W / 2.0); // floor(W/2)

    for (const auto& net : nets_hpwl) {
        std::string tag = makeNetTag(net.name);

        // ----- HardBlock pins -----
        for (int bid : net.block_ids) {
            if (bid < 0 || bid >= num_of_hardblocks) continue;
            const HardBlock& hb = hardblocks[bid];
            int i = hb.id; // assumed 0-based same as bid

            int w = get<0>(hb.block);
            int h = get<1>(hb.block);
            int dx = w / 2; // floor(w/2)
            int dy = h / 2; // floor(h/2)

            int x_idx = i + 1;
            int y_idx = i + 1;

            // max_nx >= x_i + dx  -> max_nx - x_i >= dx
            output_file << "  c" << cid++ << ": "
                        << "max_" << tag << "x - x" << x_idx
                        << " >= " << dx << "\n";

            // min_nx <= x_i + dx  -> min_nx - x_i <= dx
            output_file << "  c" << cid++ << ": "
                        << "min_" << tag << "x - x" << x_idx
                        << " <= " << dx << "\n";

            // max_ny >= y_i + dy  -> max_ny - y_i >= dy
            output_file << "  c" << cid++ << ": "
                        << "max_" << tag << "y - y" << y_idx
                        << " >= " << dy << "\n";

            // min_ny <= y_i + dy  -> min_ny - y_i <= dy
            output_file << "  c" << cid++ << ": "
                        << "min_" << tag << "y - y" << y_idx
                        << " <= " << dy << "\n";
        }

        // ----- Pad pins -----
        for (int pid : net.pad_ids) {
            // pad mapping:
            // pad0: left edge middle   -> (0, 0.5 y)
            // pad1: top edge middle    -> (floor(W/2), y)
            // pad2: right edge middle  -> (W, 0.5 y)
            // pad3: bottom edge middle -> (floor(W/2), 0)

            double padX = 0.0;

            if (pid == 0) {              // left middle
                padX = 0.0;
                // X: constant
                output_file << "  c" << cid++ << ": "
                            << "max_" << tag << "x >= " << padX << "\n";
                output_file << "  c" << cid++ << ": "
                            << "min_" << tag << "x <= " << padX << "\n";

                // Y: 0.5 y  -> max_ny >= 0.5 y; min_ny <= 0.5 y
                output_file << "  c" << cid++ << ": "
                            << "max_" << tag << "y - 0.5 y >= 0\n";
                output_file << "  c" << cid++ << ": "
                            << "min_" << tag << "y - 0.5 y <= 0\n";
            } else if (pid == 1) {       // top middle
                padX = W_mid;
                output_file << "  c" << cid++ << ": "
                            << "max_" << tag << "x >= " << padX << "\n";
                output_file << "  c" << cid++ << ": "
                            << "min_" << tag << "x <= " << padX << "\n";

                // Y: y
                output_file << "  c" << cid++ << ": "
                            << "max_" << tag << "y - y >= 0\n";
                output_file << "  c" << cid++ << ": "
                            << "min_" << tag << "y - y <= 0\n";
            } else if (pid == 2) {       // right middle
                padX = W;
                output_file << "  c" << cid++ << ": "
                            << "max_" << tag << "x >= " << padX << "\n";
                output_file << "  c" << cid++ << ": "
                            << "min_" << tag << "x <= " << padX << "\n";

                // Y: 0.5 y
                output_file << "  c" << cid++ << ": "
                            << "max_" << tag << "y - 0.5 y >= 0\n";
                output_file << "  c" << cid++ << ": "
                            << "min_" << tag << "y - 0.5 y <= 0\n";
            } else if (pid == 3) {       // bottom middle
                padX = W_mid;
                output_file << "  c" << cid++ << ": "
                            << "max_" << tag << "x >= " << padX << "\n";
                output_file << "  c" << cid++ << ": "
                            << "min_" << tag << "x <= " << padX << "\n";

                // Y: 0
                output_file << "  c" << cid++ << ": "
                            << "max_" << tag << "y >= 0\n";
                output_file << "  c" << cid++ << ": "
                            << "min_" << tag << "y <= 0\n";
            } else {
                std::cerr << "Warning: unknown pad id " << pid << " in net "
                          << net.name << "\n";
            }
        }
    }

    // ---------- Bounds / Integers / Binaries ----------
    output_file << "Bounds\n\n";
    output_file << "Integers\n\n";

    output_file << "Binaries\n  ";
    bool firstBin = true;

    // r_i binaries
    for (int i = 0; i < N; ++i) {
        if (!firstBin) output_file << " ";
        output_file << "r" << (i + 1);
        firstBin = false;
    }
    // p_ij binaries
    for (int i = 0; i < N; ++i) {
        for (int j = i + 1; j < N; ++j) {
            output_file << " p" << (i + 1) << (j + 1);
        }
    }
    // q_ij binaries
    for (int i = 0; i < N; ++i) {
        for (int j = i + 1; j < N; ++j) {
            output_file << " q" << (i + 1) << (j + 1);
        }
    }
    output_file << "\nEnd\n";

    output_file.close();
}

// -------------------------------------------------------
// main: read HB / nets / pads, call Sutanthavibul, write LP
// -------------------------------------------------------
int main(int argc, char* argv[]) {
    if (argc < 4) {
        std::cerr << "Usage: " << argv[0]
                  << " <input_file> <output_file> <x(width of die)>\n";
        return 1;
    }

    double t0, t1, t;
    t0 = GetTime(); // start time

    int num_of_hardblocks = 0;
    int num_of_pads       = 0;
    int num_of_nets       = 0;

    HardBlock* hardblocks;
    vector<vector<int>> matrix;

    std::vector<NetHPWL> nets_hpwl;

    int total_block_area = 0;

    /*---Read input---*/
    string input_filepath = argv[1];
    ifstream input_file(input_filepath);
    if (!input_file) {
        std::cerr << "Cannot open input file: " << input_filepath << "\n";
        return 1;
    }

    string lineStr;
    auto nextLogicalLine = [&](std::string& out) -> bool {
        std::string raw;
        while (std::getline(input_file, raw)) {
            auto first = raw.find_first_not_of(" \t\r\n");
            if (first == std::string::npos) continue; // all whitespace
            auto last = raw.find_last_not_of(" \t\r\n");
            out = raw.substr(first, last - first + 1);
            if (out.rfind("//", 0) == 0 || out[0] == '#') continue; // comment
            return true;
        }
        return false;
    };

    std::stringstream line;
    std::string key;

    // Read num_of_hardblocks
    if (!nextLogicalLine(lineStr)) {
        std::cerr << "Error: missing num_of_hardblocks\n";
        return 1;
    }
    line.clear(); line.str(lineStr);
    line >> key >> num_of_hardblocks;
    if (key != "num_of_hardblocks") {
        std::cerr << "Error: expected 'num_of_hardblocks', got '" << key << "'\n";
        return 1;
    }

    // Read num_of_nets
    if (!nextLogicalLine(lineStr)) {
        std::cerr << "Error: missing num_of_nets\n";
        return 1;
    }
    line.clear(); line.str(lineStr);
    line >> key >> num_of_nets;
    if (key != "num_of_nets") {
        std::cerr << "Error: expected 'num_of_nets', got '" << key << "'\n";
        return 1;
    }

    // Read num_of_pads
    if (!nextLogicalLine(lineStr)) {
        std::cerr << "Error: missing num_of_pads\n";
        return 1;
    }
    line.clear(); line.str(lineStr);
    line >> key >> num_of_pads;
    if (key != "num_of_pads") {
        std::cerr << "Error: expected 'num_of_pads', got '" << key << "'\n";
        return 1;
    }

    hardblocks = new HardBlock[num_of_hardblocks];
    total_block_area = 0;

    // Read hardblocks section: lines like "hb0 4 2" / "hbk1 413 220"
    int blocks_read = 0;
    while (blocks_read < num_of_hardblocks && nextLogicalLine(lineStr)) {
        // Skip until we see something that starts with "hb"
        // For ami33: "hbk1" 等也會通過這個檢查
        if (lineStr.rfind("hb", 0) != 0) continue;

        line.clear();
        line.str(lineStr);

        std::string name;
        int w = 0, h = 0;
        if (!(line >> name >> w >> h)) {
            std::cerr << "Error: invalid hardblock line: " << lineStr << "\n";
            return 1;
        }

        hardblocks[blocks_read].name  = name;
        hardblocks[blocks_read].block = std::make_tuple(w, h);
        hardblocks[blocks_read].id    = blocks_read;

        total_block_area += w * h;
        ++blocks_read;
    }

    if (blocks_read != num_of_hardblocks) {
        std::cerr << "Error: expected " << num_of_hardblocks
                  << " hardblocks, but read " << blocks_read << "\n";
        return 1;
    }

    // 建立「block 名稱 → index」對應，支援 hbk1/hbk5a/... 之類名字
    std::unordered_map<std::string, int> blockNameToIndex;
    for (int i = 0; i < num_of_hardblocks; ++i) {
        blockNameToIndex[hardblocks[i].name] = i;
    }

    // Read nets section: e.g.
    //   net0 hb0 hb1 hb7 pad0
    //   nDIR hbk9a hbk11
    nets_hpwl.reserve(num_of_nets);
    int nets_read = 0;
    while (nets_read < num_of_nets && nextLogicalLine(lineStr)) {
        line.clear();
        line.str(lineStr);

        NetHPWL net;
        line >> net.name;  // 可以是 net0 / nDIR / nC30 / ...

        if (net.name.empty()) continue;

        std::string tok;
        while (line >> tok) {
            // 先看是不是 hardblock 名稱（hbk1、hbk5a、hb0...）
            auto it = blockNameToIndex.find(tok);
            if (it != blockNameToIndex.end()) {
                int bid = it->second;
                net.block_ids.push_back(bid);
                continue;
            }

            // pad0~pad3
            if (tok.rfind("pad", 0) == 0) {
                int pid = std::stoi(tok.substr(3));
                if (pid < 0 || pid >= num_of_pads) {
                    std::cerr << "Error: pad index out of range in net line: "
                              << lineStr << "\n";
                    return 1;
                }
                net.pad_ids.push_back(pid);
                continue;
            }

            // 其他 token（像 GND/PWR）直接忽略，只噴 warning 一下
            if (tok != "GND" && tok != "PWR" &&
                tok != "VDD" && tok != "VSS") {
                std::cerr << "Warning: unknown token '" << tok
                          << "' in net line: " << lineStr << "\n";
            }
        }

        nets_hpwl.push_back(net);
        ++nets_read;
    }

    if (nets_read != num_of_nets) {
        std::cerr << "Error: expected " << num_of_nets
                  << " nets, but read " << nets_read << "\n";
        return 1;
    }

    // ---- Print hardblocks ----
    cout << "HardBlocks:\n";
    for (int i = 0; i < num_of_hardblocks; i++) {
        cout << "  " << hardblocks[i].name << " (hb" << i << "): "
             << get<0>(hardblocks[i].block) << " "
             << get<1>(hardblocks[i].block) << "\n";
    }
    cout << "total_block_area: " << total_block_area << "\n";

    // ---- Print nets ----
    cout << "Nets (HPWL pins):\n";
    for (const auto& net : nets_hpwl) {
        cout << "  " << net.name << " : blocks =";
        for (int bid : net.block_ids) {
            cout << " " << hardblocks[bid].name << "(hb" << bid << ")";
        }
        if (!net.pad_ids.empty()) {
            cout << " ; pads =";
            for (int pid : net.pad_ids) {
                cout << " pad" << pid;
            }
        }
        cout << "\n";
    }

    int dieWidth = std::stoi(argv[3]);

    // Build Sutanthavibul geometric constraints
    Sutanthavibul(hardblocks, num_of_hardblocks, dieWidth, matrix);

    t1 = GetTime();
    t  = t1 - t0;
    (void)t; // 如果你不想印時間，避免 unused warning

    // --- Write LP file ---
    ConvertMatrixToLPFILE(matrix, argv[2],
                          hardblocks,
                          num_of_hardblocks,
                          nets_hpwl,
                          dieWidth);

    cout << "\nWriting result to " << argv[2] << endl;

    delete [] hardblocks;

    return 0;
}
