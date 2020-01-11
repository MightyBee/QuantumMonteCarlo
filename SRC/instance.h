#include <string>
#include <vector>
#include "ConfigFile.tcc"

template int ConfigFile::get<int>(const std::string& key) const;
template double ConfigFile::get<double>(const std::string& key) const;
template unsigned int ConfigFile::get<unsigned int>(const std::string& key) const;
template size_t ConfigFile::get<size_t>(const std::string& key) const;
template std::string ConfigFile::get<std::string>(const std::string& key) const;
