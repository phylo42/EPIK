#include "exceptions.h"

#include <iostream>
#include <sstream>

//--------------------------------------------------------------------------
bad_options::bad_options(const std::string& message)
        :what_(message)
{ }

const char* bad_options::what() const throw()
{
    return std::string("Argument error" + what_).c_str();
}

//--------------------------------------------------------------------------
conflicting_options::conflicting_options(const std::string& a, const std::string& b)
        : a_(a), b_(b)
{ }

const char* conflicting_options::what() const throw()
{
    std::ostringstream message_stream;
    message_stream << "Conflicting options: " << a_ << " and " << b_ << std::endl;
    return message_stream.str().c_str();
}