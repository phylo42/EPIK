#ifndef RAPPAS_CPP_EXCEPTIONS_H
#define RAPPAS_CPP_EXCEPTIONS_H

#include <exception>
#include <string>

//--------------------------------------------------------------------------------------------------------------------
class bad_options final : public std::exception
{
public:
    bad_options(const std::string& message);
    virtual const char* what() const throw();

private:
    const std::string what_;
};

//--------------------------------------------------------------------------------------------------------------------
class conflicting_options final : public std::exception
{
public:
    conflicting_options(const std::string& a, const std::string& b);

    virtual const char* what() const throw();

private:
    const std::string a_;
    const std::string b_;
};


#endif