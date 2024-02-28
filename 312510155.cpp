#include <cstdlib>
#include <cstdio>
#include <iostream>
#include <fstream>
#include <cstring>
#include <vector>
#include <queue>
#include <map>
#include <stack>
#include <ctime>
#include <algorithm>
#include <iomanip>
using namespace std;

// enum type def for readable
enum net_t {
    IN_NET, NET, OUT_NET
};

enum gate_t {
    NOR, NAND, INV
};

enum table_index_t {
    CELL_RISE, CELL_FALL, RISE_TRANSITION, FALL_TRANSITION, RISE_POWER, FALL_POWER
};

enum pin_t : bool {
    A1, A2
};

// class predeclaration
class gate;
class wire;

// filter useless char, keep character or number or . or _
bool char_is_valid(char &c) {
    return ( ( c >= 'a' && c <= 'z' ) || ( c >= 'A' && c<= 'Z' ) || ( c >= '0' && c <= '9' ) || ( c == '.' ) || ( c == '_' ) );
}

// net class
class net{
public:
    string name;
    double transition_time;
    double load;
    double delay;
    bool value;
    bool save_value;
    gate* in_gate;
    vector<gate*> out_gate;
    net(string input_name) {
        this->name = input_name;
        this->transition_time = 0;
        this->load = 0;
        this->delay = 0;
        this->value = false;
        this->save_value = false;
        this->in_gate = nullptr;
    }
};

// gate class
class gate{
public:
    gate_t gate_type;
    string name;
    // if inverter, only A1 is used
    net* A1;
    net* A2;
    net* ZN;
    bool is_visited;
    int rise_count;
    int fall_count;
    double delay;
    double internal_power;
    gate(gate_t input_gate_type, string input_name) {
        this->gate_type = input_gate_type;
        this->name = input_name;
        this->A1 = nullptr;
        this->A2 = nullptr;
        this->ZN = nullptr;
        this->rise_count = 0;
        this->fall_count = 0;
        this->delay = 0;
        this->internal_power = 0;
        this->is_visited = false;
    }
};

// delete all useless char to get an useful word
string get_word(ifstream &input_file) {
    string word;
    char c;
    while(true)
    {
        input_file.get(c);
        // terminate condition
        if (input_file.eof()) {
            if(!word.empty()){
                return word;
            }
            else {
                word = "no_input";
                return word;
            }
        }
        // if valid char, add in word
        if(char_is_valid(c)) {
            word.push_back(c);
        }
        // else if word is not none, add to world array
        else if(!word.empty()) {
            return word;
        }
    }
}

// parse library file info
void parse_lib_info_to_array(ifstream &input_file, double (&table) [3][6][7][7], double (&cap_table) [3][2], double (&load_index) [7], double (&tran_index) [7]) {
    gate_t gate_type = NOR;
    table_index_t delay_type = CELL_RISE;
    pin_t pin_type = A1;
    bool do_delay_input = false;
    bool do_cap_input = false;
    bool do_power_input = false;
    string word;
    while(!input_file.eof()) {
        word = get_word(input_file);
        // parse gate type
        if(word == "NOR2X1") {
            gate_type = NOR;
        }
        else if(word == "INVX1") {
            gate_type = INV;
        }
        else if(word == "NANDX1") {
            gate_type = NAND;
        }
        // parse delay type
        else if(word == "cell_rise") {
            delay_type = CELL_RISE;
        }
        else if(word == "cell_fall") {
            delay_type = CELL_FALL;
        }
        else if(word == "rise_transition") {
            delay_type = RISE_TRANSITION;
        }
        else if(word == "fall_transition") {
            delay_type = FALL_TRANSITION;
        }
        // parse power type
        else if(word == "rise_power") {
            delay_type = RISE_POWER;
        }
        else if(word == "fall_power") {
            delay_type = FALL_POWER;
        }
        // parse pin type
        else if(word == "A1") {
            pin_type = A1;
            do_cap_input = true;
        }
        else if(word == "A2") {
            pin_type = A2;
            do_cap_input = true;
        }
        else if(word == "I") {
            pin_type = A1;
            do_cap_input = true;
        }
        else if(word == "ZN") {
            do_cap_input = false;
        }
        // parse timing table or power table
        else if(word == "internal_power") {
            do_delay_input = false;
            do_power_input = true;
        }
        else if(word == "timing") {
            do_delay_input = true;
            do_power_input = false;
        }
        // parse capacitance
        else if(word == "capacitance" && do_cap_input) {
            word = get_word(input_file);
            cap_table[gate_type][pin_type] = stod(word);
            if(gate_type == INV)
                cap_table[gate_type][A2] = 0;
        }
        // parse delay table
        else if(word == "values" && do_delay_input) {
            for(unsigned int i = 0; i < 7; i++) {
                for(unsigned int j = 0; j < 7; j++) {
                    word = get_word(input_file);
                    table[gate_type][delay_type][i][j] = stod(word);
                }
            }
        }
        // parse power table
        else if(word == "values" && do_power_input) {
            for(unsigned int i = 0; i < 7; i++) {
                for(unsigned int j = 0; j < 7; j++) {
                    word = get_word(input_file);
                    table[gate_type][delay_type][i][j] = stod(word);
                }
            }
        }
        // parse index
        else if(word == "index_1") {
            for(unsigned int i = 0; i < 7; i++) {
                word = get_word(input_file);
                load_index[i] = stod(word);
            }
        }
        else if(word == "index_2") {
            for(unsigned int i = 0; i < 7; i++) {
                word = get_word(input_file);
                tran_index[i] = stod(word);
            }
        }
    }
}

// split net file with keywords
bool split_net_file(string &s) {
    if(s == "output" || s == "input" || s == "wire" || s == "INVX1" || s == "NANDX1" || s == "NOR2X1" || s == "endmodule")
        return true;
    else
        return false;
}

// generate net and gate array
void parse_net_info_to_array(ifstream &input_file, ifstream &input_pat_file, vector <net*> &input_net_array, vector <gate*> &gate_array, double (&cap_table) [3][2]) {
    // split net file into different input block
    net* net_pointer;
    gate* gate_pointer;
    string word;
    bool get_new_word = true;
    // use map to find net with name as key
    map<string, net*> net_map;
    unsigned int input_net_count = 0;
    while(!input_file.eof()) {
        if(get_new_word)
            word = get_word(input_file);
        else {
            get_new_word = true;
        }
        // make inverter
        if(word == "INVX1") {
            word = get_word(input_file);
            gate_pointer = new gate(INV, word);
            gate_array.push_back(gate_pointer);
            word = get_word(input_file);
            while(!split_net_file(word)) {
                if(word == ".I") {
                    word = get_word(input_file);
                    net_pointer = net_map[word];
                    gate_pointer->A1 = net_pointer;
                    net_pointer->load += cap_table[INV][A1];
                    net_pointer->out_gate.push_back(gate_pointer);
                    word = get_word(input_file);
                }
                else if(word == ".ZN") {
                    word = get_word(input_file);
                    net_pointer = net_map[word];
                    gate_pointer->ZN = net_pointer;
                    net_pointer->in_gate = gate_pointer;
                    word = get_word(input_file);
                }
                else {
                    word = get_word(input_file);
                }
            }
            get_new_word = false;
        }
        // make nand
        else if(word == "NANDX1") {
            word = get_word(input_file);
            gate_pointer = new gate(NAND, word);
            gate_array.push_back(gate_pointer);
            word = get_word(input_file);
            while(!split_net_file(word)) {
                if(word == ".A1") {
                    word = get_word(input_file);
                    net_pointer = net_map[word];
                    gate_pointer->A1 = net_pointer;
                    net_pointer->load += cap_table[NAND][A1];
                    net_pointer->out_gate.push_back(gate_pointer);
                    word = get_word(input_file);
                }
                else if(word == ".A2") {
                    word = get_word(input_file);
                    net_pointer = net_map[word];
                    gate_pointer->A2 = net_pointer;
                    net_pointer->load += cap_table[NAND][A2];
                    net_pointer->out_gate.push_back(gate_pointer);
                    word = get_word(input_file);
                }
                else if(word == ".ZN") {
                    word = get_word(input_file);
                    net_pointer = net_map[word];
                    gate_pointer->ZN = net_pointer;
                    net_pointer->in_gate = gate_pointer;
                    word = get_word(input_file);
                }
                else {
                    word = get_word(input_file);
                }
            }
            get_new_word = false;
        }
        // make nor
        else if(word == "NOR2X1") {
            word = get_word(input_file);
            gate_pointer = new gate(NOR, word);
            gate_array.push_back(gate_pointer);
            word = get_word(input_file);
            while(!split_net_file(word)) {
                if(word == ".A1") {
                    word = get_word(input_file);
                    net_pointer = net_map[word];
                    gate_pointer->A1 = net_pointer;
                    net_pointer->load += cap_table[NOR][A1];
                    net_pointer->out_gate.push_back(gate_pointer);
                    word = get_word(input_file);
                }
                else if(word == ".A2") {
                    word = get_word(input_file);
                    net_pointer = net_map[word];
                    gate_pointer->A2 = net_pointer;
                    net_pointer->load += cap_table[NOR][A2];
                    net_pointer->out_gate.push_back(gate_pointer);
                    word = get_word(input_file);
                }
                else if(word == ".ZN") {
                    word = get_word(input_file);
                    net_pointer = net_map[word];
                    gate_pointer->ZN = net_pointer;
                    net_pointer->in_gate = gate_pointer;
                    word = get_word(input_file);
                }
                else {
                    word = get_word(input_file);
                }
            }
            get_new_word = false;
        }
        // make output net
        else if(word == "output") {
            word = get_word(input_file);
            while(!split_net_file(word)) {
                net_pointer = new net(word);
                net_pointer->load = 0.03;
                net_map[word] = net_pointer;
                word = get_word(input_file);
            }
            get_new_word = false;
        }
        // make input net
        else if(word == "input") {
            word = get_word(input_file);
            while(!split_net_file(word)) {
                net_pointer = new net(word);
                net_map[word] = net_pointer;
                input_net_count++;
                word = get_word(input_file);
            }
            get_new_word = false;
        }
        // make net
        else if(word == "wire") {
            word = get_word(input_file);
            while(!split_net_file(word)) {
                net_pointer = new net(word);
                net_map[word] = net_pointer;
                word = get_word(input_file);
            }
            get_new_word = false;
        }
    }
    while(!input_pat_file.eof() && input_net_count != 0) {
        word = get_word(input_pat_file);
        if(word != "input") {
            input_net_array.push_back(net_map[word]);
            input_net_count--;
        }
    }
}

bool parse_one_pat(ifstream &input_pat_file, vector <net*> &input_net_array) {
    string word = "x";
    unsigned int input_net_count = 0;
    while(!input_pat_file.eof() && input_net_count != input_net_array.size()) {
        word = get_word(input_pat_file);
        if(word == ".end") {
            return false;
        }
        else if(word == "1") {
            input_net_array[input_net_count]->value = true;
        }
        else {
            input_net_array[input_net_count]->value = false;
        }
        input_net_count ++;
    }
    return true;
}

// do one polation
double single_polation(double axis_0, double axis_1, double value_0, double value_1, double target_axis) {
    if(target_axis == axis_0)
        return value_0;
    else if(target_axis == axis_1)
        return value_1;
    else {
        double result = value_0 + ( (target_axis - axis_0) * (value_1 - value_0) ) / (axis_1 - axis_0);
        return result;
    }
}

// calculate delay in a while loop for queue
void calculate_delay(queue<net*> &net_queue,  double (&table) [3][6][7][7], double (&load_index) [7], double (&tran_index) [7]) {
    net* net_pointer;
    net* dominate_input_net;
    gate* gate_pointer;
    double tran;
    double cal_delay_result = 0;
    double cal_tran_result = 0;
    double cal_power_result = 0;
    while(!net_queue.empty())
    {
        net_pointer = net_queue.front();
        net_queue.pop();
        gate_pointer = net_pointer->in_gate;
        if(gate_pointer != nullptr) {
            // calculate its drive gate
            // calculate delay
            // find dominate net by value and calculate this net value
            if(gate_pointer->gate_type == INV) {
                dominate_input_net = gate_pointer->A1;
                net_pointer->value = !(dominate_input_net->value);
            }
            // nand gate
            else if(gate_pointer->gate_type == NAND) {
                // both 1
                if(gate_pointer->A1->value && gate_pointer->A2->value) {
                    net_pointer->value = false;
                    if(gate_pointer->A1->delay > gate_pointer->A2->delay) {
                        dominate_input_net = gate_pointer->A1;
                    }
                    else {
                        dominate_input_net = gate_pointer->A2;
                    }
                }
                // one 1
                else if(gate_pointer->A1->value) {
                    net_pointer->value = true;
                    dominate_input_net = gate_pointer->A2;
                }
                else if(gate_pointer->A2->value) {
                    net_pointer->value = true;
                    dominate_input_net = gate_pointer->A1;
                }
                // both 0
                else {
                    net_pointer->value = true;
                    if(gate_pointer->A1->delay > gate_pointer->A2->delay) {
                        dominate_input_net = gate_pointer->A2;
                    }
                    else {
                        dominate_input_net = gate_pointer->A1;
                    }
                }
            }
            // nor gate
            else {
                // both 1
                if(gate_pointer->A1->value && gate_pointer->A2->value) {
                    net_pointer->value = false;
                    if(gate_pointer->A1->delay < gate_pointer->A2->delay) {
                        dominate_input_net = gate_pointer->A1;
                    }
                    else {
                        dominate_input_net = gate_pointer->A2;
                    }
                }
                // one 1
                else if(gate_pointer->A1->value) {
                    net_pointer->value = false;
                    dominate_input_net = gate_pointer->A1;
                }
                else if(gate_pointer->A2->value) {
                    net_pointer->value = false;
                    dominate_input_net = gate_pointer->A2;
                }
                // both 0
                else {
                    net_pointer->value = true;
                    if(gate_pointer->A1->delay < gate_pointer->A2->delay) {
                        dominate_input_net = gate_pointer->A2;
                    }
                    else {
                        dominate_input_net = gate_pointer->A1;
                    }
                }
            }
            // find transition time
            tran = dominate_input_net->transition_time;
            // find load
            // find table index
            int c0_index;
            int c1_index;
            int t0_index;
            int t1_index;
            if(net_pointer->load < load_index[0]) {
                c0_index = 0;
                c1_index = 1;
            }
            else if(net_pointer->load > load_index[6]) {
                c0_index = 5;
                c1_index = 6;
            }
            else {
                for(c1_index = 1; c1_index < 7; c1_index++) {
                    if(net_pointer->load < load_index[c1_index])
                        break;
                }
                c0_index = c1_index - 1;
            }
            if(tran < tran_index[0]) {
                t0_index = 0;
                t1_index = 1;
            }
            else if(tran > tran_index[6]) {
                t0_index = 5;
                t1_index = 6;
            }
            else {
                for(t1_index = 1; t1_index < 7; t1_index++) {
                    if(tran < tran_index[t1_index])
                        break;
                }
                t0_index = t1_index - 1;
            }
            
            table_index_t tran_type;
            table_index_t delay_type;
            table_index_t power_type;
            if(net_pointer->value) {
                delay_type = CELL_RISE;
                tran_type = RISE_TRANSITION;
                power_type = RISE_POWER;
            }
            else {
                delay_type = CELL_FALL;
                tran_type = FALL_TRANSITION;
                power_type = FALL_POWER;
            }

            // calculate delay and transition by interpolation
            double A;
            double B;

            A = single_polation(tran_index[t0_index], tran_index[t1_index], table[gate_pointer->gate_type][delay_type][t0_index][c0_index], table[gate_pointer->gate_type][delay_type][t1_index][c0_index], tran);
            B = single_polation(tran_index[t0_index], tran_index[t1_index], table[gate_pointer->gate_type][delay_type][t0_index][c1_index], table[gate_pointer->gate_type][delay_type][t1_index][c1_index], tran);
            cal_delay_result = single_polation(load_index[c0_index], load_index[c1_index], A, B, net_pointer->load);

            A = single_polation(tran_index[t0_index], tran_index[t1_index], table[gate_pointer->gate_type][tran_type][t0_index][c0_index], table[gate_pointer->gate_type][tran_type][t1_index][c0_index], tran);
            B = single_polation(tran_index[t0_index], tran_index[t1_index], table[gate_pointer->gate_type][tran_type][t0_index][c1_index], table[gate_pointer->gate_type][tran_type][t1_index][c1_index], tran);
            cal_tran_result = single_polation(load_index[c0_index], load_index[c1_index], A, B, net_pointer->load);
            
            A = single_polation(tran_index[t0_index], tran_index[t1_index], table[gate_pointer->gate_type][power_type][t0_index][c0_index], table[gate_pointer->gate_type][power_type][t1_index][c0_index], tran);
            B = single_polation(tran_index[t0_index], tran_index[t1_index], table[gate_pointer->gate_type][power_type][t0_index][c1_index], table[gate_pointer->gate_type][power_type][t1_index][c1_index], tran);
            cal_power_result = single_polation(load_index[c0_index], load_index[c1_index], A, B, net_pointer->load);
            
            //  update gate and output net information
            net_pointer->transition_time = cal_tran_result;
            gate_pointer->delay = cal_delay_result;
            net_pointer->delay = dominate_input_net->delay + cal_delay_result;
            gate_pointer->internal_power = cal_power_result;
            // for its output gate, add them to queue if is visited berfore, cause it mean its delay can be determined
            // otherwise just reord one input delay, wait until the other delay is calculated
            for(unsigned int i = 0; i < net_pointer->out_gate.size(); i++) {
                if(net_pointer->out_gate[i]->is_visited) {
                    // is visit one time, add its out net to queue
                    net_pointer->out_gate[i]->is_visited = false;
                    net_queue.push(net_pointer->out_gate[i]->ZN);
                }
                else {
                    // not visit yet, if it is inverter, add its out net to queue
                    net_pointer->out_gate[i]->is_visited = true;
                    if(net_pointer->out_gate[i]->gate_type == INV)
                        net_queue.push(net_pointer->out_gate[i]->ZN);
                }
            }
        }
        else {
            // is input net
            for(unsigned int i = 0; i < net_pointer->out_gate.size(); i++) {
                if(net_pointer->out_gate[i]->is_visited) {
                    // is visit one time add its out net to queue
                    net_pointer->out_gate[i]->is_visited = false;
                    net_queue.push(net_pointer->out_gate[i]->ZN);
                }
                else {
                    // not visit yet, if it is inverter, add its out net to queue
                    net_pointer->out_gate[i]->is_visited = true;
                    if(net_pointer->out_gate[i]->gate_type == INV) {
                        net_queue.push(net_pointer->out_gate[i]->ZN);
                        net_pointer->out_gate[i]->is_visited = false;
                    }
                }
            }
        }
    }
}

// compare function for sort
bool compare_gate_by_name(gate* g1, gate* g2) {
    string s1 = g1->name;
    string s2 = g2->name;
    s1 = s1.erase(0,1);
    s2 = s2.erase(0,1);
    int i1 = stoi(s1);
    int i2 = stoi(s2);
    return i1 < i2;
}

// main function
int main(int argc, char** argv)
{
    ifstream net_file;
    ifstream pat_file;
    ifstream lib_file;
    
    FILE* load_output_file;
    FILE* gate_info_output_file;
    FILE* gate_power_output_file;
    FILE* coverage_output_file;

    // open input file
    net_file.open(argv[1]);
    pat_file.open(argv[2]);
    lib_file.open(argv[3]);

    // open output file
    string front_file_name = "312510155_";
    int cnt = 0;
    while(argv[1][cnt] != '.'){
        front_file_name.push_back(argv[1][cnt]);
        cnt++;
    }
    load_output_file = fopen((front_file_name + "_load.txt").c_str(), "w");
    gate_info_output_file = fopen((front_file_name + "_gate_info.txt").c_str(), "w");
    gate_power_output_file = fopen((front_file_name + "_gate_power.txt").c_str(), "w");
    coverage_output_file = fopen((front_file_name + "_coverage.txt").c_str(), "w");
    //path_output_file = fopen((front_file_name + "_path.txt").c_str(), "w");
    
    vector <string> word_array;
    // parse libery file
    double table[3][6][7][7];
    double cap_table[3][2];
    double load_index[7];
    double tran_index[7];
    parse_lib_info_to_array(lib_file, table, cap_table, load_index,tran_index);
    
    // parse net file
    vector <net*> input_net_array;
    vector <gate*> gate_array;
    parse_net_info_to_array(net_file, pat_file, input_net_array, gate_array, cap_table);
    
    // sort gate array
    sort(gate_array.begin(), gate_array.end(), compare_gate_by_name);

    // output load txt
    for(unsigned int i = 0; i < gate_array.size(); i++) {
        fprintf(load_output_file, "%s %.6f\n", gate_array[i]->name.c_str(), gate_array[i]->ZN->load);
        // calculate switch power
    }
    
    // calculate
    queue <net*> net_queue;
    int pat_count = 1;
    double total_power = 0;
    double toogle_rate = 0;
    int toogle_count = 0;
    while(parse_one_pat(pat_file, input_net_array)) {
        for(unsigned int i = 0; i < input_net_array.size(); i++) {
            net_queue.push(input_net_array[i]);
        }
        calculate_delay(net_queue,  table, load_index, tran_index);
        total_power = 0;
        toogle_count = 0;
        // calculage data for coverage txt
        for(unsigned int i = 0; i < gate_array.size(); i++) {
            // output gate info txt
            fprintf(gate_info_output_file, "%s %d %.6f %.6f\n", gate_array[i]->name.c_str(), gate_array[i]->ZN->value, gate_array[i]->delay, gate_array[i]->ZN->transition_time);
            // output gate power txt
            fprintf(gate_power_output_file, "%s %.6f %.6f\n", gate_array[i]->name.c_str(), gate_array[i]->internal_power, ((0.5) * gate_array[i]->ZN->load * 0.9 * 0.9));
            // calculate total power
            total_power += gate_array[i]->internal_power;
            if(gate_array[i]->ZN->save_value != gate_array[i]->ZN->value)
                total_power += ((0.5) * gate_array[i]->ZN->load * 0.9 * 0.9);
            // sum rise and fall toogle
            if(gate_array[i]->ZN->save_value && !gate_array[i]->ZN->value)
                gate_array[i]->fall_count += 1;
            else if(!gate_array[i]->ZN->save_value && gate_array[i]->ZN->value)
                gate_array[i]->rise_count += 1;
            // calculate coverage
            if(gate_array[i]->rise_count >= 20)
                toogle_count += 20;
            else
                toogle_count += gate_array[i]->rise_count;
            if(gate_array[i]->fall_count >= 20)
                toogle_count += 20;
            else
                toogle_count += gate_array[i]->fall_count;
            // save previous value
            gate_array[i]->ZN->save_value = gate_array[i]->ZN->value;
        }
        fprintf(gate_info_output_file, "\n");
        fprintf(gate_power_output_file, "\n");
        toogle_rate = double(100) * toogle_count / (40 * gate_array.size());
        fprintf(coverage_output_file, "%d %.6f %.2f%%\n", pat_count, total_power, toogle_rate);
        fprintf(coverage_output_file, "\n");
        pat_count++;
    }
    return 0;
}