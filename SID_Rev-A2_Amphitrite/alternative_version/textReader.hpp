#ifndef textReader_hpp
#define textReader_hpp

#include <iostream>
#include <fstream>
#include <sstream>
#include <istream>
#include <string>
#include <vector>
using namespace std;





class TextReader{

private:
	string filename;
	string foldername;


public:
	TestReader(string filename, string foldername){

		string var_list_file = "var_list.txt";
		string var_info_file = "var_info.txt";
		vector<string> var_list;
		vector<string> var_type;
		vector<string> var_unit;
		vector<float> var_val_float;
		vector<string> val_vas_str;

		var_list = textToVector(string filename, string foldername);

		for(int i = 0 , i < var_list.size(), i++){



		}

	}

	vector<string> textToVector(string filename, string foldername);

	int stringFinder(string var_name, vector<string> textfile);



};
#endif
