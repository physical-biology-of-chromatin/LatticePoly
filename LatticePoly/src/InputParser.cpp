//
//  InputParser.cpp
//  LatticePoly
//
//  Created by mtortora on 21/12/2019.
//  Adapted from https://www.dreamincode.net/forums/topic/183191-create-a-simple-configuration-file-parser/

//  Copyright © 2019 ENS Lyon. All rights reserved.
//

#include <fstream>
#include <sstream>
#include <typeinfo>

#include "InputParser.hpp"


int Nrelax;
int Nmeas;
int Ninter;

int Nchain;

int NliqMC;
int Ndrop;

bool InitDrop;
bool RestartFromFile;

double Kint;

double R;
double Ldens;
double Rconfinement;

double Jll;
double Jlp;
double Jlpp;
double Jpp;
double Jppp;
double Jns;
double Jpppp;
double EV;

double inactiveRatio;
double propRate;

int propagationMode;

double painterAct;
double cisSpread;
double transSpread;
double readerWriter;
double boost;
double nucleoTurn;

double originRate;
double replicRate;

std::string latticeType;
std::string polyType;

std::string outputDir;
std::string domainPath;
std::string painterPath;


InputParser::InputParser(const std::string& _filePath): filePath(_filePath)
{
	ExtractKeys();
}

void InputParser::ParseVars()
{
	Nrelax          = GetValueOfKey<int>("Nrelax");
	Nmeas           = GetValueOfKey<int>("Nmeas");
	Ninter          = GetValueOfKey<int>("Ninter");
	
	Nchain          = GetValueOfKey<int>("Nchain");

	NliqMC          = GetValueOfKey<int>("NliqMC");
	Ndrop           = GetValueOfKey<int>("Ndrop");
	
	InitDrop        = GetValueOfKey<bool>("InitDrop");
	RestartFromFile = GetValueOfKey<bool>("RestartFromFile");

	Kint            = GetValueOfKey<double>("Kint");
	
	R               = GetValueOfKey<double>("R");
	Ldens           = GetValueOfKey<double>("Ldens");
	Rconfinement    = GetValueOfKey<double>("Rconfinement");

	Jll             = GetValueOfKey<double>("Jll");
	Jlp             = GetValueOfKey<double>("Jlp");
    Jlpp            = GetValueOfKey<double>("Jlpp");
	Jpp             = GetValueOfKey<double>("Jpp");
    Jns             = GetValueOfKey<double>("Jns");
    Jppp            = GetValueOfKey<double>("Jppp");
	Jpppp           = GetValueOfKey<double>("Jpppp");
	EV              = GetValueOfKey<double>("EV");

	inactiveRatio   = GetValueOfKey<double>("inactiveRatio");
	propRate        = GetValueOfKey<double>("propRate");

	propagationMode = GetValueOfKey<int>("propagationMode");

	painterAct      = GetValueOfKey<double>("painterAct");
	cisSpread       = GetValueOfKey<double>("cisSpread");
	transSpread     = GetValueOfKey<double>("transSpread");
	readerWriter    = GetValueOfKey<double>("readerWriter");
	boost           = GetValueOfKey<double>("boost");
	nucleoTurn      = GetValueOfKey<double>("nucleoTurn");

	originRate      = GetValueOfKey<double>("originRate");
	replicRate      = GetValueOfKey<double>("replicRate");
	
	polyType        = GetValueOfKey<std::string>("polyType");
	latticeType     = GetValueOfKey<std::string>("latticeType");
	
	outputDir       = GetValueOfKey<std::string>("outputDir");
	domainPath      = GetValueOfKey<std::string>("domainPath");
	painterPath     = GetValueOfKey<std::string>("painterPath");
}

void InputParser::ExtractKeys()
{
	std::ifstream file(filePath);
		
	if ( !file.good() )
		throw std::runtime_error("InputParser: Couldn't open input file " + filePath);

	std::string line;
	
	size_t lineNo = 0;
	
	while ( std::getline(file, line) )
	{
		std::string tmp = line;

		if ( tmp.empty() )
			continue;

		RemoveComment(tmp);
		
		if ( OnlyWhitespace(tmp) )
			continue;

		ParseLine(tmp, ++lineNo);
	}

	file.close();
}

void InputParser::ExtractContents(const std::string& line)
{
	std::string key;
	std::string value;
	
	std::string tmp = line;
	
	tmp.erase(0, tmp.find_first_not_of("\t "));
	
	size_t sepPos = tmp.find('=');
	
	ExtractKey(key, sepPos, tmp);
	ExtractValue(value, sepPos, tmp);

	if ( !KeyExists(key) )
		contents.insert(std::pair<std::string, std::string>(key, value));
	else
		throw std::runtime_error("InputParser: Can only have unique key names");
}

void InputParser::ExtractKey(std::string& key, size_t const& sepPos, const std::string& line) const
{
	key = line.substr(0, sepPos);
	
	if ( key.find('\t') != line.npos || key.find(' ') != line.npos )
		key.erase(key.find_first_of("\t "));
}

void InputParser::ExtractValue(std::string& value, size_t const& sepPos, const std::string& line) const
{
	value = line.substr(sepPos + 1);
	
	value.erase(0, value.find_first_not_of("\t "));
	value.erase(value.find_last_not_of("\t ") + 1);
}

void InputParser::RemoveComment(std::string& line) const
{
	if ( line.find(';') != line.npos )
		line.erase(line.find(';'));
}

void InputParser::ParseLine(const std::string& line, size_t const lineNo)
{
	if ( line.find('=') == line.npos )
		throw std::runtime_error("InputParser: Couldn't find separator on line: " + Converter::T_to_string(lineNo));
	
	if ( !ValidLine(line) )
		throw std::runtime_error("InputParser: Bad format for line: " + Converter::T_to_string(lineNo));

	ExtractContents(line);
}

bool InputParser::ValidLine(const std::string& line) const
{
	std::string tmp = line;
	
	tmp.erase(0, tmp.find_first_not_of("\t "));
	
	if ( tmp[0] == '=' )
		return false;

	for ( size_t i = tmp.find('=') + 1; i < tmp.length(); ++i )
	{
		if ( tmp[i] != ' ' )
			return true;
	}

	return false;
}

bool InputParser::OnlyWhitespace(const std::string& line) const
{
	return ( line.find_first_not_of(' ') == line.npos );
}

bool InputParser::KeyExists(const std::string& key) const
{
	return contents.find(key) != contents.end();
}

template <typename ValueType>
ValueType InputParser::GetValueOfKey(const std::string& key) const
{
	if ( !KeyExists(key) )
		throw std::runtime_error("No entry found for input parameter " + key + " in file " + filePath);
	
	return Converter::string_to_T<ValueType>(contents.find(key)->second);
}

template<typename T>
std::string InputParser::Converter::T_to_string(T const& val)
{
	std::ostringstream ostr;
	ostr << val;

	return ostr.str();
}

template<typename T>
T InputParser::Converter::string_to_T(std::string const& val)
{
	std::istringstream istr(val);
	
	T returnVal;
	
	if ( !(istr >> returnVal) )
		throw std::runtime_error("InputParser: Not a valid " + (std::string) typeid(T).name());

	return returnVal;
}
