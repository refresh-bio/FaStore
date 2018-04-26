/*
  This file is a part of FaStore software distributed under GNU GPL 2 licence.
  The code in this file is based on ORCOM software.

  Github:	https://github.com/refresh-bio/FaStore

  Authors: Lukasz Roguski, Idoia Ochoa, Mikel Hernaez & Sebastian Deorowicz
*/

#ifndef H_EXCEPTION
#define H_EXCEPTION

#include "Globals.h"

#include <string>
#include <stdexcept>


class Exception : public std::exception
{
	std::string message;

public:
	Exception(const char* msg_)
		: message(msg_)
	{}

	Exception(const std::string& msg_)
		: message(msg_)
	{}

	~Exception() throw()
	{}

	const char* what() const throw()				// for std::exception interface
	{
		return message.c_str();
	}
};


#endif // H_EXCEPTION
