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
int NG1;

int Nchain;

int NliqMC;
int Ndrop;

bool InitDrop;
bool RestartFromFile;
bool neigh;

double Kint;

double R;
double Ldens;

double Jll;
double Jlp;
double Jpp;

double Jaa;
double Jbb;
double Jab;

double Jtad_a;
double Jtad_b;



double Jf;
double Jf_sister;
double Jpair;
int Centromere;
double stop_replication;

double inactiveRatio;
double propRate;
double keco1;

int propagationMode;
int cohesionMode;
int ForkTableMode;
double originRate;
double replicRate;
int Ndf;
int enhancement_fork;
int enhancement_sister;
int enhancement_cohesin;

double N_extruders;
double permeability;
double loading_rate;
double unloading_rate;
bool StartFromPODLS;
double n_barriers;

std::string latticeType;
std::string polyType;

std::string outputDir;
std::string domainPath;
std::string colorPath;
std::string CARpath;
std::string PODLSpath;
std::string OriginsPath;





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
	NG1          = GetValueOfKey<int>("NG1");
	Centromere          = GetValueOfKey<int>("Centromere");

	
	NliqMC          = GetValueOfKey<int>("NliqMC");
	Ndrop           = GetValueOfKey<int>("Ndrop");
	
	InitDrop        = GetValueOfKey<bool>("InitDrop");
	RestartFromFile = GetValueOfKey<bool>("RestartFromFile");
	StartFromPODLS = GetValueOfKey<bool>("StartFromPODLS");


	
	Kint            = GetValueOfKey<double>("Kint");
	
	R               = GetValueOfKey<double>("R");
	Ldens           = GetValueOfKey<double>("Ldens");
	stop_replication           = GetValueOfKey<double>("stop_replication");

	
	Jll             = GetValueOfKey<double>("Jll");
	Jlp             = GetValueOfKey<double>("Jlp");
	Jpp             = GetValueOfKey<double>("Jpp");
	Jpair             = GetValueOfKey<double>("Jpair");
	
	
	
	Jaa             = GetValueOfKey<double>("Jaa");
	Jbb             = GetValueOfKey<double>("Jbb");
	Jab             = GetValueOfKey<double>("Jab");

	Jtad_a             = GetValueOfKey<double>("Jtad_a");
	Jtad_b             = GetValueOfKey<double>("Jtad_b");
	
	Jf              = GetValueOfKey<double>("Jf");
	Jf_sister              = GetValueOfKey<double>("Jf_sister");

	Ndf              = GetValueOfKey<int>("Ndf");
	neigh = GetValueOfKey<bool>("neigh");
	enhancement_fork = GetValueOfKey<int>("enhancement_fork");
	enhancement_sister = GetValueOfKey<int>("enhancement_sister");
	enhancement_cohesin = GetValueOfKey<int>("enhancement_cohesin");


	N_extruders        = GetValueOfKey<double>("N_extruders");
	permeability		=GetValueOfKey<double>("permeability");
	loading_rate		=GetValueOfKey<double>("loading_rate");
	unloading_rate		=GetValueOfKey<double>("unloading_rate");
	n_barriers          = GetValueOfKey<double>("n_barriers");




	
	inactiveRatio   = GetValueOfKey<double>("inactiveRatio");
	keco1   = GetValueOfKey<double>("keco1");

	propRate        = GetValueOfKey<double>("propRate");

	propagationMode = GetValueOfKey<int>("propagationMode");
	cohesionMode =  GetValueOfKey<int>("cohesionMode");
	ForkTableMode =  GetValueOfKey<int>("ForkTableMode");

	originRate      = GetValueOfKey<double>("originRate");
	replicRate      = GetValueOfKey<double>("replicRate");
	
	polyType        = GetValueOfKey<std::string>("polyType");
	latticeType     = GetValueOfKey<std::string>("latticeType");
	
	outputDir       = GetValueOfKey<std::string>("outputDir");
	
	domainPath      = GetValueOfKey<std::string>("domainPath");
	colorPath     = GetValueOfKey<std::string>("colorPath");
	CARpath     = GetValueOfKey<std::string>("CARpath");
	PODLSpath     = GetValueOfKey<std::string>("PODLSpath");
	OriginsPath		=GetValueOfKey<std::string>("OriginsPath");


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
