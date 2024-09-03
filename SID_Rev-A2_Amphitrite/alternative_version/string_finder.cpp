#include <iostream>
#include <fstream>
#include <sstream>
#include <istream>
#include <string>
#include <vector>
using namespace std;

vector<string> textToVector(string filename){

  vector<string> textfile;

  ifstream infile("test_file.txt");
  string line;

  while (std::getline(infile, line)){
    istringstream iss(line);
    textfile.push_back(line);
  }
  return textfile;
}



int stringFinder(string var_name, vector<string> textfile){

  int lineNumber = textfile.size();
  int max_it = 20;

  for(int i = 0; i < lineNumber; i++){

    if(textfile[i] == var_name){

      return i;

    }
  }
}

int main(){

  string filename = "test_file.txt";
  string var_name = "test_name";
  vector<string> textfile;
  int var_line;
  float var_value;

  textfile = textToVector(filename);

  var_line = stringFinder(var_name, textfile);

  cout << "Searched var is in line " << var_line+1 << " (var: " << textfile[var_line] << ")" << endl;

  var_value = stof(textfile[var_line+1]);

  cout << var_name << " = " << var_value << endl;

}
