//
// Created by meital on 12/17/19.
//

#ifndef COMPARISON_SMALLFIELDS_H
#define COMPARISON_SMALLFIELDS_H


#include "NTL/ZZ_p.h"
#include "NTL/ZZ.h"
#ifdef __x86_64__
#include <x86intrin.h>
#elif __aarch64__
#include "../infra/SSE2NEON.h"
#endif
#include <gmp.h>
#include <libscapi/include/primitives/Prg.hpp>
#include <libscapi/include/primitives/Mersenne.hpp>


using namespace std;
using namespace NTL;

class ZpMersenne13 {


//private:
public: //TODO return to private after tesing

    static const unsigned int p = 8191;
    unsigned int elem;

public:

    ZpMersenne13(){elem = 0;};
    ZpMersenne13(long elem){ this->elem = elem %p;}//no need to optimize, just an int modulo



    ZpMersenne13& operator=(const ZpMersenne13& other){elem = other.elem; return *this;};
    bool operator!=(const ZpMersenne13& other){ return !(other.elem == elem); };
    bool operator==(const ZpMersenne13& other){ return other.elem == elem; };

    ZpMersenne13 operator+(const ZpMersenne13& f2)
    {
        ZpMersenne13 answer;

        answer.elem = (elem + f2.elem);

        if(answer.elem>=p)
            answer.elem-=p;

        return answer;
    }
    ZpMersenne13 operator-(const ZpMersenne13& f2)
    {
        ZpMersenne13 answer;

        int temp =  (int)elem - (int)f2.elem;

        if(temp<0){
            answer.elem = temp + p;
        }
        else{
            answer.elem = temp;
        }

        return answer;
    }
    ZpMersenne13 operator/(const ZpMersenne13& f2)
    {
        //code taken from NTL for the function XGCD
        int a = f2.elem;
        int b = p;
        long s;

        int  u, v, q, r;
        long u0, v0, u1, v1, u2, v2;

        int aneg = 0, bneg = 0;

        if (a < 0) {
            if (a < -NTL_MAX_LONG) Error("XGCD: integer overflow");
            a = -a;
            aneg = 1;
        }

        if (b < 0) {
            if (b < -NTL_MAX_LONG) Error("XGCD: integer overflow");
            b = -b;
            bneg = 1;
        }

        u1=1; v1=0;
        u2=0; v2=1;
        u = a; v = b;

        while (v != 0) {
            q = u / v;
            r = u % v;
            u = v;
            v = r;
            u0 = u2;
            v0 = v2;
            u2 =  u1 - q*u2;
            v2 = v1- q*v2;
            u1 = u0;
            v1 = v0;
        }

        if (aneg)
            u1 = -u1;


        s = u1;

        if (s < 0)
            s =  s + p;

        ZpMersenne13 inverse(s);

        return inverse* (*this);
    }

    ZpMersenne13 operator*(const ZpMersenne13& f2)
    {
        ZpMersenne13 answer;


        //the multiplication can be put into an int, since it has no more than 26 bits.
        answer.elem = (elem * f2.elem)%p;

        return answer;
    }

    ZpMersenne13& operator+=(const ZpMersenne13& f2){
        elem = (f2.elem + elem);

        if(elem>=p)
            elem-=p;

        return *this;
    };
    ZpMersenne13& operator*=(const ZpMersenne13& f2)
    {
        elem = (elem * f2.elem) %p;

        return *this;
    }

    ZpMersenne13 sqrt()
    {
        //The algorithm for checking the square root of a value is as follows:
        //We know that 2^31 and 2^61 are both divisible by 4 (the results are 2^29 and 2^59 respectively). So 2^31-1=3 mod 4 and 2^61-1=3 mod 4.
        //So if we have b=x^2 (over Mersenne61) then we can compute x by b^{2^59}.
        //To do this, we can make about 58 field multiplications:
        //Set b_1 = b, then
        //For i=2...59:
        //compute b_i = (b_{i-1})^2.
        //So x1=b_59 and x2=-b_59 = 2^61-1-b_59
        //Check that x1^2 = b, if it does then output it, otherwise, it means that a cheat is detected.
        ZpMersenne13 answer = *this;
        for (int i=2; i<=12; i++){
            answer *= answer;
        }
        ZpMersenne13 check = answer*answer;

        if (check != *this){
            cout<<"CHEATING!!!"<<endl;
            return ZpMersenne13(0);
        }

        return answer;
    }

};

inline ::ostream& operator<<(::ostream& s, const ZpMersenne13& a){ return s << a.elem; };

class Zp16BitPrime {


//private:
public: //TODO return to private after tesing

    static const unsigned int p = 65521;
    unsigned int elem;

public:

    Zp16BitPrime(){elem = 0;};
    Zp16BitPrime(long elem){ this->elem = elem %p;}//no need to optimize, just an int modulo



    Zp16BitPrime& operator=(const Zp16BitPrime& other){elem = other.elem; return *this;};
    bool operator!=(const Zp16BitPrime& other){ return !(other.elem == elem); };
    bool operator==(const Zp16BitPrime& other){ return other.elem == elem; };

    Zp16BitPrime operator+(const Zp16BitPrime& f2)
    {
        Zp16BitPrime answer;

        answer.elem = (elem + f2.elem);

        if(answer.elem>=p)
            answer.elem-=p;

        return answer;
    }
    Zp16BitPrime operator-(const Zp16BitPrime& f2)
    {
        Zp16BitPrime answer;

        int temp =  (int)elem - (int)f2.elem;

        if(temp<0){
            answer.elem = temp + p;
        }
        else{
            answer.elem = temp;
        }

        return answer;
    }
    Zp16BitPrime operator/(const Zp16BitPrime& f2)
    {
        //code taken from NTL for the function XGCD
        int a = f2.elem;
        int b = p;
        long s;

        int  u, v, q, r;
        long u0, v0, u1, v1, u2, v2;

        int aneg = 0, bneg = 0;

        if (a < 0) {
            if (a < -NTL_MAX_LONG) Error("XGCD: integer overflow");
            a = -a;
            aneg = 1;
        }

        if (b < 0) {
            if (b < -NTL_MAX_LONG) Error("XGCD: integer overflow");
            b = -b;
            bneg = 1;
        }

        u1=1; v1=0;
        u2=0; v2=1;
        u = a; v = b;

        while (v != 0) {
            q = u / v;
            r = u % v;
            u = v;
            v = r;
            u0 = u2;
            v0 = v2;
            u2 =  u1 - q*u2;
            v2 = v1- q*v2;
            u1 = u0;
            v1 = v0;
        }

        if (aneg)
            u1 = -u1;


        s = u1;

        if (s < 0)
            s =  s + p;

        Zp16BitPrime inverse(s);

        return inverse* (*this);
    }

    Zp16BitPrime operator*(const Zp16BitPrime& f2)
    {
        Zp16BitPrime answer;


        //the multiplication can be put into an int, since it has no more than 26 bits.
        answer.elem = (elem * f2.elem)%p;

        return answer;
    }

    Zp16BitPrime& operator+=(const Zp16BitPrime& f2){
        elem = (f2.elem + elem);

        if(elem>=p)
            elem-=p;

        return *this;
    };
    Zp16BitPrime& operator*=(const Zp16BitPrime& f2)
    {
        elem = (elem * f2.elem) %p;

        return *this;
    }

    Zp16BitPrime sqrt()
    {
        //The algorithm for checking the square root of a value is as follows:
        //We know that 2^31 and 2^61 are both divisible by 4 (the results are 2^29 and 2^59 respectively). So 2^31-1=3 mod 4 and 2^61-1=3 mod 4.
        //So if we have b=x^2 (over Mersenne61) then we can compute x by b^{2^59}.
        //To do this, we can make about 58 field multiplications:
        //Set b_1 = b, then
        //For i=2...59:
        //compute b_i = (b_{i-1})^2.
        //So x1=b_59 and x2=-b_59 = 2^61-1-b_59
        //Check that x1^2 = b, if it does then output it, otherwise, it means that a cheat is detected.
        Zp16BitPrime answer = *this;
        for (int i=2; i<=12; i++){
            answer *= answer;
        }
        Zp16BitPrime check = answer*answer;

        if (check != *this){
            cout<<"CHEATING!!!"<<endl;
            return Zp16BitPrime(0);
        }

        return answer;
    }

};

inline ::ostream& operator<<(::ostream& s, const Zp16BitPrime& a){ return s << a.elem; };



#endif //COMPARISON_SMALLFIELDS_H
