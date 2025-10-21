// Inlcudes
#include <fstream>
#include <iostream>
#include <string>
#include <vector>
#include <limits>
#include <stdexcept>
#include <cctype>
#include <cmath>
#include <algorithm>

int findNumberOfCollums(std::fstream &file, char delim = ';') {

    std::string line = "";
    if(!std::getline(file, line)) return 0;
    if(line.empty()) return 0;

    int collums = 1;
    for(char c: line) {
        if(c == '\r') continue;
        if(c == delim) collums++;
    }

    file.clear();
    file.seekg(0);
    return collums;
}

int findNumberOfRows(std::fstream &file) {
    int rows = 0;
    std::string empty;
    
    while(std::getline(file, empty)) rows++;

    file.clear();
    file.seekg(0);
    return rows;
}

std::vector<std::string>splitOnDelim(const std::string &s, const int nCollums ,char delim = ';') {
    std::vector<std::string>cells; std::string current;
    cells.reserve(nCollums);

    for(char c : s) {
        if(c == delim) {
            cells.push_back(current);
            current.clear();
        } else {
            current.push_back(c);
        }
    }
    cells.push_back(current);

    return cells;
}

std::vector<std::vector<int>>calcSpaces(const std::vector<std::vector<std::string>>table,const int nCollums) {
    std::vector<std::string>header = table[0];
    std::vector<std::vector<int>>result;

    for(std::vector<std::string>row : table) {
        std::vector<int>t_Coll;
        t_Coll.reserve(nCollums);
        for(size_t i = 0; i < nCollums; i++) {
            int currentLen = row[i].length();
            int res = header[i].length() - currentLen;

            t_Coll.push_back(res);
        }
        result.push_back(t_Coll);
    }

    return result;
}

// Main Function
int main(int argc, char* argv[]) {
    std::string path = argv[argc - 1]; 
    std::fstream File(path);

    if(!File.is_open()) return -1;
    
    int nCollums = findNumberOfCollums(File);
    int nRows = findNumberOfRows(File);

    std::cout << "Collums: " << nCollums << std::endl;
    std::cout << "Rows: " << nRows << std::endl;

    std::string line;
    std::vector<std::vector<std::string>>rows;
    rows.reserve(nRows);

    File.clear();
    File.seekg(0);
    while(std::getline(File, line)) rows.push_back(splitOnDelim(line, nCollums));
    File.close();

    if(rows.size() == 0) return -1;
    
    std::cout << "Content:" << std::endl;
    
    std::vector<std::vector<int>>nWhitespaces = calcSpaces(rows, nCollums);
    for(size_t r = 0; r < rows.size(); r++) {
        std::cout << '|';
        for(size_t c  = 0; c < rows.data()->size(); c++) {
            std::cout << rows[r][c];
            for(int i = 0; i < nWhitespaces[r][c]; i++) std::cout << " "; 
            std::cout << '|';
        }
        std::cout << std::endl; 
    }
    
    return 0;
}
