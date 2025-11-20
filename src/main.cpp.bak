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
#include "datastructure.hpp"
#define EQ 0
#define GE 1
#define LE 2
using namespace std;

// output: matrix, core util

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

// input: hardblocks[0..N-1] with (w_i, h_i)
// output: matrix with 2N^2 rows and ( (N+1)^2 + 1 ) columns:
//         variables (x,y,r,p,q, global y), OP, RHS
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

    auto idx_p = [N, numPermu](int i, int j) {
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
            //     -> y_i - y_j + r_j(h_j - w_j) - M p_ij - M q_ij >= -2M + h_j
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

void ConvertMatrixToLPFILE(vector<vector<int>>& matrix, string output_filepath) {
    ofstream output_file(output_filepath);
    if (!output_file) return;
    if (matrix.empty()) { output_file.close(); return; }

    int numRows = (int)matrix.size();
    int numCols = (int)matrix[0].size();

    int RHScol = numCols - 1;
    int OPcol  = numCols - 2;
    int ycol   = numCols - 3; // global y

    // deduce N from #constraints: 2N^2 = numRows
    int N = (int)std::lround(std::sqrt(numRows / 2.0));
    int numPermu = N * (N - 1) / 2;

    auto idx_x = [N](int i) { return i; };
    auto idx_y = [N](int i) { return N + i; };
    auto idx_r = [N](int i) { return 2 * N + i; };

    auto idx_p = [N, numPermu](int i, int j) {
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
    output_file << "  y\n";

    // ---------- Subject To ----------
    output_file << "Subject To\n";

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

int main(int argc, char* argv[]) {
    if (argc < 4) {
        std::cerr << "Usage: " << argv[0] << " <input_file> <output_file> <x(width of die)\n";
        return 1;
    }

    double t0, t1, t;
    t0 = GetTime(); // start time

    int num_of_hardblocks; HardBlock* hardblocks;
    int num_of_pads;       Pad* pads;
    int num_of_nets;       Net* nets;
    vector<vector<int>> matrix;
    int total_block_area = 0;
 

/*---Read input---*/
    string input_filepath = argv[1];  // ../testcase/public1.txt
    ifstream input_file(input_filepath);
    string input;
    stringstream line;

    // Read Number of Hardblocks
    getline(input_file, input);
    string dummy0;
    line.clear(); line.str(input);
    line >> dummy0 >> num_of_hardblocks;
    hardblocks = new HardBlock[num_of_hardblocks];

    // 不要再手動多讀幾行，直接由 while 來跳過空白 / 註解
    total_block_area = 0;

    int blocks_read = 0;
    while (blocks_read < num_of_hardblocks && std::getline(input_file, input)) {
        // 去掉前後空白
        auto first = input.find_first_not_of(" \t\r\n");
        if (first == std::string::npos) continue; // 全空白
        auto last = input.find_last_not_of(" \t\r\n");
        std::string trimmed = input.substr(first, last - first + 1);

        // 跳過註解
        if (trimmed.rfind("//", 0) == 0 || trimmed[0] == '#') continue;

        line.clear();
        line.str(trimmed);

        std::string name;
        int w = 0, h = 0;
        if (!(line >> name >> w >> h)) {
            std::cerr << "Error: invalid hardblock line: " << trimmed << "\n";
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

        
    // pirnt hardblocks
    cout << "HardBlocks: " << endl;
    for (int i = 0; i < num_of_hardblocks; i++) {
        cout << "HardBlock " << i << ": " << get<0>(hardblocks[i].block) << " " << get<1>(hardblocks[i].block) << endl;
    }
    cout << "total_block_area: " << total_block_area << endl;


    Sutanthavibul(hardblocks, num_of_hardblocks, std::stoi(argv[3]), matrix);


    t1 = GetTime();
    t = t1 - t0;

// /* --- Write result to output file --- */
    ConvertMatrixToLPFILE(matrix, argv[2]);
    cout << endl << "Writing result to " << argv[2] << endl;
 
//     ofstream output_file(output_filepath);
//     output_file << "Wirelength " << total_wirelength << endl << endl;
//     output_file << "NumHardBlocks " << num_of_hardblocks << endl;
//     for (int i = 0; i < num_of_hardblocks; i++) {
//         output_file << "hb" << i << " " << get<0>(cg.coord[i]) << " " << get<1>(cg.coord[i]) << " " << cg.rotate[i] << endl;
//     }
//     output_file.close();
// 
//     double Stage2_time = GetTime() - t0 - Stage1_time;
//     cout << "Stage 1 time: " << Stage1_time << " seconds" << endl;
//     cout << "Stage 2 time: " << Stage2_time << " seconds" << endl;
//     cout << "Total   Time: " << GetTime() - t0 << " seconds" << endl;

    return 0;
} 
