//
// Created by meital on 12/17/19.
//

#include "SmallFields.h"


template <>
TemplateField<ZpMersenne13>::TemplateField(long fieldParam) {

    this->fieldParam = 8191;
    this->elementSizeInBytes = 2;//round up to int. no need to optimize space and communication here
    this->elementSizeInBits = 13;

    auto randomKey = prg.generateKey(128);
    prg.setKey(randomKey);

    m_ZERO = new ZpMersenne13(0);
    m_ONE = new ZpMersenne13(1);
}

template <>
ZpMersenne13 TemplateField<ZpMersenne13>::GetElement(long b) {


    if(b == 1)
    {
        return *m_ONE;
    }
    if(b == 0)
    {
        return *m_ZERO;
    }
    else{
        ZpMersenne13 element(b);
        return element;
    }
}

template <>
void TemplateField<ZpMersenne13>::elementToBytes(unsigned char* elemenetInBytes, ZpMersenne13& element){

    memcpy(elemenetInBytes, (byte*)(&element.elem), 2);
}


template <>
void TemplateField<ZpMersenne13>::elementVectorToByteVector(vector<ZpMersenne13> &elementVector, vector<byte> &byteVector){

    copy_byte_array_to_byte_vector((byte *)elementVector.data(), elementVector.size()*elementSizeInBytes, byteVector,0);
}


template <>
ZpMersenne13 TemplateField<ZpMersenne13>::bytesToElement(unsigned char* elemenetInBytes){

    return ZpMersenne13((unsigned short int)(*(unsigned short int *)elemenetInBytes));
}


template <>
TemplateField<Zp16BitPrime>::TemplateField(long fieldParam) {

    this->fieldParam = 65521;
    this->elementSizeInBytes = 4;//round up to int. no need to optimize space and communication here
    this->elementSizeInBits = 16;

    auto randomKey = prg.generateKey(128);
    prg.setKey(randomKey);

    m_ZERO = new Zp16BitPrime(0);
    m_ONE = new Zp16BitPrime(1);
}

template <>
Zp16BitPrime TemplateField<Zp16BitPrime>::GetElement(long b) {


    if(b == 1)
    {
        return *m_ONE;
    }
    if(b == 0)
    {
        return *m_ZERO;
    }
    else{
        Zp16BitPrime element(b);
        return element;
    }
}

template <>
void TemplateField<Zp16BitPrime>::elementToBytes(unsigned char* elemenetInBytes, Zp16BitPrime& element){

    memcpy(elemenetInBytes, (byte*)(&element.elem), 4);
}


template <>
void TemplateField<Zp16BitPrime>::elementVectorToByteVector(vector<Zp16BitPrime> &elementVector, vector<byte> &byteVector){

    copy_byte_array_to_byte_vector((byte *)elementVector.data(), elementVector.size()*elementSizeInBytes, byteVector,0);
}


template <>
Zp16BitPrime TemplateField<Zp16BitPrime>::bytesToElement(unsigned char* elemenetInBytes){

    return Zp16BitPrime((unsigned int)(*(unsigned int *)elemenetInBytes));
}