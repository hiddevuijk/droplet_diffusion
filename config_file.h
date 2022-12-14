
#ifndef GUARD_CONFIG_FILE_H
#define GUARD_CONFIG_FILE_H

#include <string>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <map>
#include <sstream>

class Config {
public:
	Config( const char * fname,
			const char separator = '=',
			const char comment = '#');


	template<class T>
	T get_parameter(const std::string& parameterName) const
	{
		// check if parameter name is in dict
		T val;
		std::stringstream ss(dict.find(parameterName)->second);
		ss >> val;
		
		return val;	

	}
private:
	std::map< std::string, std::string > dict;


};

namespace configFunc {
	void trim(std::string& str);
};



Config::Config(const char * fname, const char separator_char ,
               const char comment_char )
{

	// read file
	std::ifstream file(fname);

	// check if file exists
	if ( !file.is_open() )
		throw "fuck this ";

	std::string line, lhs, rhs;
	size_t cmnt_index, sep_index;

	//for(; getline(file,line); ) {
	while( std::getline(file,line) ) {
		

		// check for comment, if # is inline, remove rhs
		cmnt_index = line.find(comment_char);
		if(cmnt_index != std::string::npos) {
			line.erase(cmnt_index, std::string::npos);
		}
		
		if( line.size() == 0) continue;

		// fine the separator char
		sep_index = line.find(separator_char);
	
		// if not found skip this line
		if( sep_index != std::string::npos ) {
			// check if valid	
			lhs = line.substr(0,sep_index);
			rhs = line.substr(sep_index+1,line.size()-1);
			//remove leading and trailing white spaces
			configFunc::trim(lhs);
			configFunc::trim(rhs);		
			// check if not empty!!!
			// check for double keys!!!
			dict[lhs] = rhs;	
		}
	}


}

void configFunc::trim(std::string& str)
{

	unsigned int N = str.size();
	// check if N > 0

	// remove leading white space 
	while(true) {
		if( str[0] == ' ') {
			str.erase(0,1);
		} else {
			break;
		}
	}

	// remove trailing white space
	for( unsigned int i=N-1;i>0;--i) {
		if(str[i] == '\0') continue;
		if(str[i] == ' ' ) {
			str.erase(str.end()-1);
		} else {
			break;
		}
	}
}


#endif // GUARD_CONFIGFILE_H
