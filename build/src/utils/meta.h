#ifndef RAPPAS_CPP_TEMPLATES_H
#define RAPPAS_CPP_TEMPLATES_H

template <bool flag, class IsTrue, class IsFalse>
struct choose;

template <class IsTrue, class IsFalse>
struct choose<false, IsTrue, IsFalse> {
    typedef IsFalse type;
};

template <class IsTrue, class IsFalse>
struct choose<true, IsTrue, IsFalse> {
    typedef IsTrue type;
};

#endif //RAPPAS_CPP_TEMPLATES_H
