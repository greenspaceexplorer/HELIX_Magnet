#ifndef TOKEN_PARSER_H
#define TOKEN_PARSER_H

#include <iterator>
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>

using namespace std;

class token_parser{
    public:
        token_parser(char token_str = ','){ token = token_str; }
        string const &operator[](size_t index) const{
            return row_vec[index];
        }
        size_t size(){ return row_vec.size(); }
        void read_row(istream &csv_file){
            string line;
            getline(csv_file,line);

            stringstream stream_line(line);
            string cell;

            row_vec.clear();
            
            while(getline(stream_line,cell,token)){
                row_vec.push_back(cell);
            }
            if (!stream_line && cell.empty()){
                row_vec.push_back("");
            }
        }
    private:
        vector<string> row_vec;
        char token;
};

istream &operator>>(istream &str,token_parser& row){
       row.read_row(str);
       return str;
};

#endif // TOKEN_PARSER_H
