//
//  InputParser.hpp
//  LatticePoly
//
//  Created by mtortora on 21/12/2019.
//  Copyright © 2019 ENS Lyon. All rights reserved.
//

#ifndef InputParser_hpp
#define InputParser_hpp

#include <map>

#include "globals.hpp"


class InputParser
{
public:
	InputParser(const std::string&);
	
	void ParseVars();
	
private:
	std::string filePath;
	std::map<std::string, std::string> contents;

	void RemoveComment(std::string&) const;
	void ParseLine(const std::string&, size_t const);
	
	void ExtractKey(std::string&, size_t const&, const std::string&) const;
	void ExtractValue(std::string&, size_t const&, const std::string&) const;
	void ExtractContents(const std::string&);

	void ExtractKeys();

	bool OnlyWhitespace(const std::string&) const;
	bool ValidLine(const std::string&) const;
	bool KeyExists(const std::string&) const;

	template <typename ValueType>
	ValueType GetValueOfKey(const std::string&) const;
	
	struct Converter
	{
		template<typename T>
		static std::string T_to_string(T const&);
		
		template<typename T>
		static T string_to_T(std::string const&);
	};
};


#endif /* InputParser_hpp */
