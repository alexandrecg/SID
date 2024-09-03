#include <iostream>
#include <fstream>
#include <sstream>
#include <istream>
#include <string>
#include <vector>

using namespace std;

//Extract First Float
float eff(string str)
{
    stringstream ss;

    /* Storing the whole string into string stream */
    ss << str;

    /* Running loop till the end of the stream */
    string temp;
    float found;
    while (!ss.eof()) {

        /* extracting word by word from stream */
        ss >> temp;

        /* Checking the given word is integer or not */
        if (stringstream(temp) >> found)
            return found;

        /* To save from space at the end of string */
        temp = "";
    }
}

int main(){

  int c = 0;
  int max_it = 20;

  vector<float> read_values;

  ifstream infile("test_file.txt");
  string line;

  while (std::getline(infile, line) or c > max_it){
      istringstream iss(line);

      read_values.push_back(eff(line));

      cout << "Reading line: " << line << endl;
      cout << "Reading value: " << eff(line) << endl;
      cout << "Counter value: " << c << endl;
      c++;
  }

  cout << "\nvector size: " << read_values.size() << endl;
  cout << "\nReading the entire vector: " << endl;


  for(int i = 0; i < read_values.size(); i++){

    cout << read_values[i] << endl;

  }

}
