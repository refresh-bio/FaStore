/*
  This file is a part of FaStore software distributed under GNU GPL 2 licence.

  Github:	https://github.com/refresh-bio/FaStore

  Authors: Lukasz Roguski, Idoia Ochoa, Mikel Hernaez & Sebastian Deorowicz
*/

#ifndef VERSION_H
#define VERSION_H

#include "fastore_bin/Globals.h"
#include <string>

// WARN: version.cpp file should be automated automatically!!!

std::string GetAppVersion();
std::string GetCompilationTime();

#endif // VERSION_H
