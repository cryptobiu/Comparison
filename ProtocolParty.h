//
// Created by moriya on 12/5/19.
//

#ifndef COMPARISON_PROTOCOLPARTY_H
#define COMPARISON_PROTOCOLPARTY_H

#include <vector>
#include <libscapi/include/primitives/Mersenne.hpp>
#include <libscapi/include/cryptoInfra/Protocol.hpp>
#include <libscapi/include/primitives/Prg.hpp>
#include <libscapi/include/primitives/Matrix.hpp>

using namespace std;

#define flag_print false
#define flag_print_timings true

template <class FieldType>
class ProtocolParty : public MPCProtocol, public HonestMajority {

private:
    int N, T, partyId;//number of servers

    int times; //number of times to run the run function
    int iteration; //number of the current iteration

    int lambda, delta;
    vector<FieldType> twoSquares;

    boost::asio::io_service io_service;
    vector<PrgFromOpenSSLAES> prgs;

    TemplateField<FieldType> *field;
    vector<shared_ptr<ProtocolPartyData>> parties;
    Measurement* timer;

    vector<FieldType> alpha; // N distinct non-zero field elements
    vector<FieldType> beta;
    VDM<FieldType> matrix_vand;
    HIM<FieldType> matrix_him;
    VDMTranspose<FieldType> matrix_vand_transpose;
    HIM<FieldType> matrix_for_interpolate;
    vector<FieldType> y_for_interpolate;

    vector<vector<FieldType>> polynomials;
    vector<FieldType> randomTAnd2TShares;
    int randomTAnd2TSharesOffset = 0;
    vector<FieldType> randomSharesArray;
    int randomSharesOffset = 0;

    vector<FieldType> randomsToLSB;
    vector<vector<FieldType>> bitsToLSB;
    int randomWithBitsOffset;


    void initializationPhase();
    bool preparationPhase();
    void interpolatePolynomial(vector<FieldType>& coeff, const FieldType* a, FieldType* b, int size);
    FieldType evalPolynomial(const vector<FieldType>& coeff, vector<FieldType> & a);

    bool RandomSharingForInputs(int no_random, vector<FieldType> & sharingBufInputsTElements);
    void generateRandom2TAndTShares(int numOfRandomPairs, vector<FieldType>& randomElementsToFill);
    void generateRandomShares(int numOfRandoms, vector<FieldType> &randomElementsToFill);

    void calcSendBufElements(vector<vector<FieldType>> & sendBufsElements, PrgFromOpenSSLAES & prg, int start, int end);
    void calcRecBufElements(vector<vector<FieldType>> & recBufsElements, vector<FieldType> & randomElementsToFill, int start, int end);

    void openShare(vector<FieldType> &shares, vector<FieldType> &secrets, int d);

    void openShareSetRecBuf(vector<FieldType> &shares, vector<FieldType> &secrets, int d, vector<vector<byte>> &recBufsBytes);

    void roundFunctionSync(vector<vector<byte>> &sendBufs, vector<vector<byte>> &recBufs, int round);
    void exchangeData(vector<vector<byte>> &sendBufs,vector<vector<byte>> &recBufs, int first, int last);
    void roundFunctionSyncElements(vector<vector<FieldType>> &sendBufs, vector<vector<FieldType>> &recBufs, int round);
    void exchangeDataElements(vector<vector<FieldType>> &sendBufs,vector<vector<FieldType>> &recBufs, int first, int last);

    /**
     * The cheap way: Create a HIM from the αi’s onto ZERO (this is actually a row vector), and multiply
     * this HIM with the given x-vector (this is actually a scalar product).
     * The first (and only) element of the output vector is the secret.
     */
    FieldType interpolate(vector<FieldType>& x);
    FieldType reconstructShare(vector<FieldType>& x, int d);

    void DNHonestMultiplication(FieldType *a, FieldType *b, vector<FieldType> &cToFill, int numOfTupples);
    void inverse(FieldType *a, vector<FieldType> &bToFill);

    void randomBit(vector<FieldType> &bToFill);
    void checkRandomBit(vector<FieldType> &bits);

    FieldType unboundedOR(FieldType* shares, int l);
    void checkUnboundedOR(FieldType* bits, int numBits);

    void prefixOR(FieldType* shares, vector<FieldType> & bToFill);

    FieldType bitwiseLessThan(FieldType* a, FieldType* b, int size);

    void geneateRandomWithBits(int size, vector<FieldType> & randomToFill, vector<vector<FieldType>> & bitsToFill, int numRandoms);
    void getBits(FieldType x, vector<FieldType> & bits, int size);
    FieldType LSB(FieldType & a, int bitsSize);

    FieldType compare(FieldType & a, FieldType & b, int bitSize);

    void comparePhase();

public:
    ProtocolParty(int argc, char* argv[]);

    /**
     * This method runs the protocol:
     * 1. Preparation Phase
     * 2. Input Phase
     * 3. Computation Phase
     * 4. Verification Phase
     * 5. Output Phase
     */
    void run() override;

    bool hasOffline() override {
        return true;
    }


    bool hasOnline() override {
        return true;
    }

    /**
     * This method runs the protocol:
     * Preparation Phase
     */
    void runOffline() override;

    /**
     * This method runs the protocol:
     * Input Phase
     * Computation Phase
     * Verification Phase
     * Output Phase
     */
    void runOnline() override;

};


template <class FieldType>
ProtocolParty<FieldType>::ProtocolParty(int argc, char* argv[]) : MPCProtocol("Comparison", argc, argv, false)
{

    partyId = stoi(this->getParser().getValueByKey(arguments, "partyID"));

    numParties = stoi(this->getParser().getValueByKey(arguments, "partiesNumber"));
    numThreads = stoi(this->getParser().getValueByKey(arguments, "numThreads"));

    string fieldType = this->getParser().getValueByKey(arguments, "fieldType");

    this->times = stoi(this->getParser().getValueByKey(arguments, "internalIterationsNumber"));

    vector<string> subTaskNames{"Offline", "preparationPhase", "Online", "inputPhase", "ComputePhase", "outputPhase"};
    timer = new Measurement(*this, subTaskNames);

    if(fieldType.compare("ZpMersenne31") == 0) {
        field = new TemplateField<FieldType>(2147483647);
    } else if(fieldType.compare("ZpMersenne61") == 0) {
        field = new TemplateField<FieldType>(0);
    }

    N = numParties;
    T = (numParties+1)/2 - 1;

    delta = 8;
    lambda = 4;

    MPCCommunication comm;
    string partiesFile = this->getParser().getValueByKey(arguments, "partiesFile");
    parties = comm.setCommunication(io_service, partyId, N, partiesFile);


    prgs.resize(numThreads);
    int* keyBytes = new int[4];
    for (int i=0; i<numThreads; i++){
        for (int j=0; j<4; j++){
            keyBytes[j] = field->Random().elem;
        }
        SecretKey key((byte*)keyBytes, 16, "");
        prgs[i].setKey(key);
    }
    delete [] keyBytes;

    auto t1 = high_resolution_clock::now();
    initializationPhase();

    auto t2 = high_resolution_clock::now();

    auto duration = duration_cast<milliseconds>(t2-t1).count();
    if(flag_print_timings) {
        cout << "time in milliseconds initializationPhase: " << duration << endl;
    }

}


template <class FieldType>
void ProtocolParty<FieldType>::run() {

    for (iteration=0; iteration<times; iteration++){

        auto t1start = high_resolution_clock::now();
        timer->startSubTask("Offline", iteration);
        runOffline();
        timer->endSubTask("Offline", iteration);
        timer->startSubTask("Online", iteration);
        runOnline();
        timer->endSubTask("Online", iteration);

        auto t2end = high_resolution_clock::now();
        auto duration = duration_cast<milliseconds>(t2end-t1start).count();

        cout << "time in milliseconds for protocol: " << duration << endl;
    }


}

template <class FieldType>
void ProtocolParty<FieldType>::runOffline() {
    auto t1 = high_resolution_clock::now();
    timer->startSubTask("preparationPhase", iteration);
    if(preparationPhase() == false) {
        if(flag_print) {
            cout << "cheating!!!" << '\n';}
        return;
    }
    else {
        if(flag_print) {
            cout << "no cheating!!!" << '\n' << "finish Preparation Phase" << '\n';}
    }
    timer->endSubTask("preparationPhase", iteration);
    auto t2 = high_resolution_clock::now();

    auto duration = duration_cast<milliseconds>(t2-t1).count();
    if(flag_print_timings) {
        cout << "time in milliseconds preparationPhase: " << duration << endl;
    }
}

template <class FieldType>
void ProtocolParty<FieldType>::runOnline() {



    auto t1 = high_resolution_clock::now();
    timer->startSubTask("inputPhase", iteration);
    comparePhase();
    timer->endSubTask("inputPhase", iteration);
    auto t2 = high_resolution_clock::now();
    auto duration = duration_cast<milliseconds>(t2-t1).count();

    if(flag_print_timings) {
        cout << "time in milliseconds inputPhase: " << duration << endl;
    }

//    cout<<"*********** check inverse *************"<<endl;
//    vector<FieldType> element(10);
//    vector<FieldType> inverseArr(10);
//    for (int i=0; i<10; i++) {
//        element[i] = randomSharesArray[randomSharesOffset++];
//    }
//
//    cout<<"share of random elements: "<<endl;
//    for (int i=0; i<10; i++) {
//        cout<<element[i]<<" ";
//    }
//    cout<<endl;
//
//    vector<FieldType> openedElement(10);
//    openShare(element, openedElement, T);
//    cout<<"opened elements: "<<endl;
//    for (int i=0; i<10; i++) {
//        cout<<openedElement[i]<<" ";
//    }
//    cout<<endl;
//
//    inverse(element.data(), inverseArr);
//    cout<<"elements inverse: " << endl;
//    for (int i=0; i<10; i++) {
//        cout<<inverseArr[i]<<" ";
//    }
//    cout<<endl;
//
//    vector<FieldType> openedInverse(10);
//    openShare(inverseArr, openedInverse, T);
//    cout<<"opened elements inverse: "<<endl;
//    for (int i=0; i<10; i++) {
//        cout<<openedInverse[i]<<" ";
//    }
//    cout<<endl;
//
//
//    for (int i=0; i<10; i++) {
//        auto res = openedElement[i] * openedInverse[i];
//
//        cout << "the output is " << res << endl;
//    }
//
//    cout<<"*********** check RANDOM BIT *************"<<endl;
//
//    vector<FieldType> bits(8);
//    checkRandomBit(bits);
//
//    cout<<"*********** check UNBOUNDED FAN-IN OR*************"<<endl;
//
//    checkUnboundedOR(bits.data(), 8);
//    checkUnboundedOR(bits.data(), 2);
//    checkUnboundedOR(bits.data(), 4);
//    checkUnboundedOR(bits.data(), 6);
//
//    cout<<"*********** check PREFIX OR*************"<<endl;
//
//    vector<FieldType> allBits(31);
//    vector<FieldType> allBitsPrefix(31);
//    vector<FieldType> allBitsPrefixOpened(31);
//    randomBit(allBits);
//
//    openShare(allBits, allBitsPrefixOpened, T);
//
//    cout<<"all bits:"<<endl;
//    for (int i=0; i<31; i++){
//        cout<<allBitsPrefixOpened[i]<<" ";
//    }
//    cout<<endl;
//
//    prefixOR(allBits.data(), allBitsPrefix);
//    openShare(allBitsPrefix, allBitsPrefixOpened, T);
//
//    cout<<"prefix bits:"<<endl;
//    for (int i=0; i<31; i++){
//        cout<<allBitsPrefixOpened[i]<<" ";
//    }
//    cout<<endl;
//
//    cout<<"*********** check BITWISE LESS THAN*************"<<endl;
//    vector<FieldType> a(31);
//    vector<FieldType> b(31);
//    vector<FieldType> aOpened(31);
//    vector<FieldType> bOpened(31);
//    randomBit(a);
//    randomBit(b);
//
//    openShare(a, aOpened, T);
//    cout<<"a:"<<endl;
//    for (int i=0; i<31; i++){
//        cout<<aOpened[i]<<" ";
//    }
//    cout<<endl;
//
//    openShare(b, bOpened, T);
//    cout<<"b:"<<endl;
//    for (int i=0; i<31; i++){
//        cout<<bOpened[i]<<" ";
//    }
//    cout<<endl;
//
//    vector<FieldType> resV(1);
//    resV[0] = bitwiseLessThan(a.data(), bOpened.data(), a.size());
//
//    vector<FieldType> resVOpened(1);
//    openShare(resV, resVOpened, T);
//
//    cout<<"a < b ? "<<resVOpened[0]<<endl;
//
//    cout<<"*********** check GENERATE RANDOM WITH BITS *************"<<endl;
//
//    int size = 31;
//    vector<vector<FieldType>> bitsOfRand;
//    vector<FieldType> rand(1);
//    vector<FieldType> randOpened(1);
//    geneateRandomWithBits(size, rand, bitsOfRand, 1);
//
//    vector<FieldType> bitsOfRandOpened(size);
//    openShare(bitsOfRand[0], bitsOfRandOpened, T);
//    cout<<"bits of random:"<<endl;
//    for (int i=0; i<size; i++){
//        cout<<bitsOfRandOpened[i];
//    }
//    cout<<endl;
//
//    openShare(rand, randOpened, T);
//    cout<<"actual random = "<<randOpened[0]<<endl;
//
//    cout<<"*********** check LSB *************"<<endl;
//
//    vector<FieldType> lsb(1);
//    vector<FieldType> lsbOpened(1);
//    lsb[0] = LSB(rand[0], size);
//
//    openShare(lsb, lsbOpened, T);
//    cout<<"lsb of random = "<<lsbOpened[0]<<endl;


//    cout<<"*********** check COMPARE *************"<<endl;
//
//    FieldType aC, bC, res;
//    int numCompares = 100;
//    long totalTime = 0;
//    for (int i=0; i<numCompares; i++) {
//        auto start = high_resolution_clock::now();
//        aC = randomSharesArray[randomSharesOffset++];
//        bC = randomSharesArray[randomSharesOffset++];
//
//        res = compare(aC, bC, 31);
//        auto end = high_resolution_clock::now();
//        auto duration = duration_cast<milliseconds>(end- start).count();
//        cout << "compare took " << duration << " ms"<<endl;
//        totalTime += duration;
//    }
//
//
//    vector<FieldType> temp(3);
//    temp[0] = aC;
//    temp[1] = bC;
//    temp[2] = res;
//    openShare(temp, temp, T);
//    cout<<"a = "<<temp[0]<<endl;
//    cout<<"b = "<<temp[1]<<endl;
//    cout<<"a < b ? = "<<temp[2]<<endl;
//
//    cout << "compute " <<numCompares<<" compares took " << totalTime/numCompares << " ms in average"<<endl;


}

template <class FieldType>
void ProtocolParty<FieldType>::comparePhase(){

    cout<<"*********** check COMPARE *************"<<endl;

    FieldType aC, bC, res;
    int numCompares = 100;
    long totalTime = 0;
    for (int i=0; i<numCompares; i++) {
        auto start = high_resolution_clock::now();
        aC = randomSharesArray[randomSharesOffset++];
        bC = randomSharesArray[randomSharesOffset++];

        res = compare(aC, bC, 31);
        auto end = high_resolution_clock::now();
        auto duration = duration_cast<milliseconds>(end- start).count();
//        cout << "compare took " << duration << " ms"<<endl;
        totalTime += duration;
    }


    vector<FieldType> temp(3);
    temp[0] = aC;
    temp[1] = bC;
    temp[2] = res;
    openShare(temp, temp, T);
    cout<<"a = "<<temp[0]<<endl;
    cout<<"b = "<<temp[1]<<endl;
    cout<<"a < b ? = "<<temp[2]<<endl;

    cout << "compute " <<numCompares<<" compares took " << totalTime/numCompares << " ms in average"<<endl;

}

template <class FieldType>
void ProtocolParty<FieldType>::checkRandomBit(vector<FieldType> & bits){
    int numBits = bits.size();
    vector<FieldType> openBits(numBits);
    randomBit(bits);

    cout<<"shares of random bits:"<<endl;
    for (int i=0; i<numBits; i++){
        cout<<bits[i]<<" ";
    }
    cout<<endl;
    openShare(bits, openBits, T);
    cout<<"actual random bits:"<<endl;
    for (int i=0; i<numBits; i++){
        cout<<openBits[i]<<" ";
    }
    cout<<endl;
}

template <class FieldType>
void ProtocolParty<FieldType>::checkUnboundedOR(FieldType* bits, int numBits){
    vector<FieldType> orVal(1);
    vector<FieldType> orValOpened(1);
    orVal[0] = unboundedOR(bits, numBits);
    openShare(orVal, orValOpened, T);

    cout<<"result of unbounded fan-in or of "<<numBits<<" bits = "<<orValOpened[0]<<endl;
}


template <class FieldType>
void ProtocolParty<FieldType>::inverse(FieldType *a, vector<FieldType> &bToFill) {

    int numOfElements = bToFill.size();
    vector<FieldType> raShares(numOfElements);//hold both in the same vector to send in one batch
    vector<FieldType> openedRAShares(numOfElements);
    FieldType one(1);


//    DNHonestMultiplication(a, randomSharesArray.data() + randomSharesOffset, raShares, numOfElements);
    //generate the shares for x+a and y+b. do it in the same array to send once
    for (int k = 0; k < numOfElements; k++)//go over only the logit gates
    {
        //compute the share of r*a
        raShares[k] = a[k]*randomSharesArray[randomSharesOffset + k];
//        cout<<"r*a share = "<<raShares[k]<<endl;
    }

    openShare(raShares, openedRAShares, 2*T);
//    cout<<"ra open = "<<openedRAShares[0]<<endl;
    for (int k = 0; k < numOfElements; k++)//go over only the logit gates
    {
        //compute the share of r*a
//        cout<<"inverse of r = "<< (one / openedRAShares[k])<<endl;
        bToFill[k] = randomSharesArray[randomSharesOffset + k] / openedRAShares[k];
//        cout<<"inverse of a share = "<< bToFill[k]<<endl;
    }

    randomSharesOffset += numOfElements;

}


template <class FieldType>
void ProtocolParty<FieldType>::randomBit(vector<FieldType> &bToFill) {

    int numOfElements = bToFill.size();
//    cout<<"numOfElements = "<<numOfElements<<endl;
    vector<FieldType> r(numOfElements);
    vector<FieldType> rSquareShares(numOfElements);//hold both in the same vector to send in one batch
    vector<FieldType> openedRSquareShares(numOfElements);
    FieldType one(1);
    FieldType twoInv(2);
    twoInv = one / twoInv;

//    cout<<"randomSharesOffset = "<<randomSharesOffset<<endl;
    //generate the shares for x+a and y+b. do it in the same array to send once

//    DNHonestMultiplication(randomSharesArray.data() + randomSharesOffset, randomSharesArray.data() + randomSharesOffset, rSquareShares, numOfElements);
    for (int k = 0; k < numOfElements; k++)//go over only the logit gates
    {
        //compute the share of r*a
        r[k] = randomSharesArray[randomSharesOffset + k];
        rSquareShares[k] = r[k]*r[k];
    }

//    cout<<"share r:"<<endl;
//    for (int i=0; i<numOfElements; i++){
//        cout<<r[i]<<" ";
//    }
//    cout<<endl;
//    DNHonestMultiplication(r.data(), r.data(), rSquareShares, numOfElements);

//    cout<<"share r^2:"<<endl;
//    for (int i=0; i<numOfElements; i++){
//        cout<<rSquareShares[i]<<" ";
//    }
//    cout<<endl;

    openShare(rSquareShares, openedRSquareShares, 2*T);
//    cout<<"opened r^2:"<<endl;
//    for (int i=0; i<numOfElements; i++){
//        cout<<openedRSquareShares[i]<<" ";
//    }
//    cout<<endl;

    FieldType temp;
    vector<FieldType> roots(numOfElements);
    for (int k = 0; k < numOfElements; k++)//go over only the logit gates
    {
        roots[k] = openedRSquareShares[k].sqrt();
        temp = *(field->GetZero()) - roots[k];

        if (temp.elem < roots[k].elem){
            roots[k] = temp;
        }

    }

//    cout<<"roots:"<<endl;
//    for (int i=0; i<numOfElements; i++){
//        cout<<roots[i]<<" ";
//    }
//    cout<<endl;

//    vector<FieldType> rInverse(numOfElements);
//    inverse(r.data(), rInverse);
//    cout<<"inverse r:"<<endl;
//    for (int i=0; i<numOfElements; i++){
//        cout<<rInverse[i]<<" ";
//    }
//    cout<<endl;
    for (int k = 0; k < numOfElements; k++)
    {
        bToFill[k] = r[k] / roots[k];
        bToFill[k] = (bToFill[k] + one) * twoInv;
    }
//    cout<<"bToFill first:"<<endl;
//    for (int i=0; i<numOfElements; i++){
//        cout<<bToFill[i]<<" ";
//    }
//    cout<<endl;
//
//
//    for (int k = 0; k < numOfElements; k++)
//    {
////        bToFill[k] = roots[k] * rInverse[k];
//        bToFill[k] = (bToFill[k] + one) * twoInv;
//    }
//
//    cout<<"bToFill second:"<<endl;
//    for (int i=0; i<numOfElements; i++){
//        cout<<bToFill[i]<<" ";
//    }
//    cout<<endl;

    randomSharesOffset += numOfElements;

//    openShare(rInverse, rInverse, T);
//    cout<<"opened inverse:"<<endl;
//    for (int i=0; i<numOfElements; i++){
//        cout<<rInverse[i]<<" ";
//    }
//    cout<<endl;

//    openShare(r, r, T);
//    cout<<"opened r:"<<endl;
//    for (int i=0; i<numOfElements; i++){
//        cout<<r[i]<<" ";
//    }
//    cout<<endl;

}


template <class FieldType>
FieldType ProtocolParty<FieldType>::unboundedOR(FieldType* shares, int l) {
    vector<FieldType> A(l);
    vector<FieldType> b(l);

    A[0] = *field->GetOne();
    for (int i=0; i<l; i++){
        A[0] += shares[i];
    }

    b[0] = *field->GetOne();
    for (int i=1; i<l; i++){
        A[i] = A[0];
        b[i] = randomSharesArray[randomSharesOffset + i-1];
    }

    FieldType lastB = randomSharesArray[randomSharesOffset + l-1];
    int originalOffset = randomSharesOffset;
    randomSharesOffset += l;

    //interpolate done in initialization phase
    vector<FieldType> randomsInverse(l);
    inverse(randomSharesArray.data() + originalOffset, randomsInverse);

    vector<FieldType> c(l);
    //calc first column - A * all B values
    DNHonestMultiplication(A.data(), b.data(), c, l);
    c[0] = A[0]; //the first value should be fixed

    //calc second column - A * all B values
    DNHonestMultiplication(c.data(), randomsInverse.data(), c, l);

    //open C
    vector<FieldType> cOpened(l);
    openShare(c, cOpened, T);

//    cout<<"opened c:"<<endl;
//    for (int i=0; i<l; i++){
//        cout<<cOpened[i]<<" ";
//    }
//    cout<<endl;

    FieldType cMult = cOpened[0];
    //compute all A squares
    for (int i=1; i<l-1; i++){
        cMult *= cOpened[i];
        A[i] = cMult * b[i+1];
    }
    A[l-1] = cMult * cOpened[l-1] * lastB;

//    openShare(A, cOpened, T);
//    cout<<"opened A:"<<endl;
//    for (int i=0; i<l; i++){
//        cout<<cOpened[i]<<" ";
//    }
//    cout<<endl;



    return evalPolynomial(polynomials[l-2], A);
}

template <class FieldType>
void ProtocolParty<FieldType>::prefixOR(FieldType* shares, vector<FieldType> & bToFill) {
    int size = bToFill.size();
//    cout<<"original size = "<<size<<endl;


    //calculate xi = unbounded or of each delta bits
    vector<FieldType> x(lambda);
    vector<FieldType> y(lambda);
    vector<FieldType> f(lambda);
    vector<FieldType> s(lambda);

    FieldType* sharesAlligned;
    int newSize;
    int remain = size % delta;
    if (remain == 0){
//        sharesAlligned.resize(size);
        newSize = size;
        sharesAlligned = shares;
//        memcpy((byte*)sharesAlligned.data(), (byte*)shares, size*field->getElementSizeInBytes());
    } else {
        newSize = (size/delta + 1)*delta;
        sharesAlligned = new FieldType[newSize];
//        sharesAlligned.resize((size/delta + 1)*delta);
//        cout<<"sharesAlligned size = "<<newSize<<endl;

        memcpy((byte*)sharesAlligned, (byte*)shares, size*field->getElementSizeInBytes());
        memset((byte*)sharesAlligned + size*field->getElementSizeInBytes(), 0, (newSize - size)*field->getElementSizeInBytes());
    }

//    vector<FieldType> temp(sharesAlligned.size());
//    openShare(sharesAlligned,temp,T);
//    cout<<"sharesAlligned:"<<endl;
//    for (int i=0; i<sharesAlligned.size(); i++){
//        cout<<temp[i]<<" ";
//    }
//    cout<<endl;

    for (int i=0; i<lambda; i++){
        x[i] = unboundedOR(sharesAlligned + i*delta, delta);
    }

    //calculate yi = prefix or of x
    y[0] = x[0];
    for (int i=1; i<lambda; i++){
        y[i] = unboundedOR(x.data(), i+1);
    }

    //calculate fi = find the group with the first "1" bit
    f[0] = y[0];
    for (int i=1; i<lambda; i++){
        f[i] = y[i] - y[i-1];
    }

    //find the real b values for the group which contains the first "1" bit.
    //Note that the group index is unknown
    vector<FieldType> a(delta, *field->GetZero());
    vector<FieldType> aToMultiply(lambda);

    //TODO can be optimized!!
    for (int j=0; j<delta; j++){
        for (int i=0; i<lambda; i++) {
            aToMultiply[i] = sharesAlligned[i * delta + j];
        }

        DNHonestMultiplication(f.data(), aToMultiply.data(), aToMultiply, lambda);
        for (int i=0; i<lambda; i++) {
            a[j] += aToMultiply[i];
        }
    }

    //calculate the real b values for the group which contains the first "1" bit.
    //Note that the group index is unknown
    vector<FieldType> b(delta);
    b[0] = a[0];
    for (int i=1; i<delta; i++){
        b[i] = unboundedOR(a.data(), i+1);
    }

    //calculate si = all zeros including the group with the first "1" bit and then all ones.
    for (int i=0; i<lambda; i++){
        s[i] = y[i] - f[i];

    }

    vector<FieldType> fToMultiply(delta);
    //calculate the real b values:
    for (int i=0; i<lambda; i++){
        for (int j=0; j<delta; j++){
            fToMultiply[j] = f[i];
        }

        DNHonestMultiplication(fToMultiply.data(), b.data(), fToMultiply, delta);

//        cout<<"fToMultiply:"<<endl;
//        for (int i=0; i<delta; i++){
//            cout<<fToMultiply[i]<<" ";
//        }
//        cout<<endl;
        if (i<lambda-1) {
            for (int j = 0; j < delta; j++) {
                bToFill[i * delta + j] = fToMultiply[j] + s[i];
            }
        } else { //last group
            for (int j = 0; j < remain; j++) {
                bToFill[i * delta + j] = fToMultiply[j] + s[i];
            }
        }
    }

//    openShare(x,x,T);
//    cout<<"x:"<<endl;
//    for (int i=0; i<lambda; i++){
//        cout<<x[i]<<" ";
//    }
//    cout<<endl;
//    openShare(y,y,T);
//    cout<<"y:"<<endl;
//    for (int i=0; i<lambda; i++){
//        cout<<y[i]<<" ";
//    }
//    cout<<endl;
//    openShare(f,f,T);
//    cout<<"f:"<<endl;
//    for (int i=0; i<lambda; i++){
//        cout<<f[i]<<" ";
//    }
//    cout<<endl;
//
//    openShare(b,b,T);
//    cout<<"b_i0:"<<endl;
//    for (int i=0; i<delta; i++){
//        cout<<b[i]<<" ";
//    }
//    cout<<endl;
//
//    openShare(s,s,T);
//    cout<<"s:"<<endl;
//    for (int i=0; i<lambda; i++){
//        cout<<s[i]<<" ";
//    }
//    cout<<endl;


}


template <class FieldType>
FieldType ProtocolParty<FieldType>::bitwiseLessThan(FieldType* a, FieldType* b, int size) {
//    vector<FieldType> mults(size);
//    DNHonestMultiplication(a, b, mults, size);
    FieldType two(2);

    //calculate [ci] = [ai]^[bi]
    vector<FieldType> c(size);
    for (int i=0; i<size; i++){
        c[i] = a[i] + b[i];
        c[i] = c[i] - (two * a[i] * b[i]);
    }

//    vector<FieldType> temp(size);
//    openShare(c,temp,T);
//    cout<<"c:"<<endl;
//    for (int i=0; i<size; i++){
//        cout<<temp[i]<<" ";
//    }
//    cout<<endl;

    //compute d[i] = prefix OR of ci
    vector<FieldType> d(size);
    prefixOR(c.data(), d);
//    openShare(d,temp,T);
//    cout<<"d:"<<endl;
//    for (int i=0; i<size; i++){
//        cout<<temp[i]<<" ";
//    }
//    cout<<endl;

    //compute e = [di]-[di+1]
    vector<FieldType> e(size);
    e[0] = d[0];
    for (int i=1; i<size; i++){
        e[i] = d[i] - d[i-1];
    }
//    openShare(e,temp,T);
//    cout<<"e:"<<endl;
//    for (int i=0; i<size; i++){
//        cout<<temp[i]<<" ";
//    }
//    cout<<endl;

    FieldType res(0);
    vector<FieldType> mults(size);
    DNHonestMultiplication(e.data(), b, mults, size);
    for (int i=0; i<size; i++){
//        res += e[i] * b[i];
        res += mults[i];
    }

    return res;
}

template <class FieldType>
void ProtocolParty<FieldType>::geneateRandomWithBits(int size, vector<FieldType> & randomToFill, vector<vector<FieldType>> & bitsToFill, int numRandoms) {

    randomToFill.resize(numRandoms);
    bitsToFill.resize(numRandoms);

    for (int j=0; j<numRandoms; j++) {
        bitsToFill[j].resize(size);
        randomBit(bitsToFill[j]);

        randomToFill[j] = *field->GetZero();
        for (int i = 0; i < size; i++) {
            randomToFill[j] += bitsToFill[j][i] * twoSquares[twoSquares.size() - size + i];
        }
    }
}

template <class FieldType>
FieldType ProtocolParty<FieldType>::LSB(FieldType & x, int bitsSize) {

    vector<FieldType> rBits = bitsToLSB[randomWithBitsOffset];
    FieldType r = randomsToLSB[randomWithBitsOffset++];

//    geneateRandomWithBits(bitsSize, r, rBits);


//    vector<FieldType> temp(bitsSize);
//    openShare(rBits, temp, T);
//    cout<<"bits of random:"<<endl;
//    for (int i=0; i<bitsSize; i++){
//        cout<<temp[i];
//    }
//    cout<<endl;
//    temp[0] = r;
//    openShare(temp, temp, T);
//    cout<<"r:"<<temp[0]<<endl;


    vector<FieldType> c(1);
    vector<FieldType> cOpened(1);
    c[0] = x + r;
    openShare(c, cOpened, T);
//    cout<<"c:"<<cOpened[0]<<endl;

    vector<FieldType> cBits(bitsSize);
    getBits(cOpened[0], cBits, bitsSize);


//    openShare(cBits, temp, T);
//    cout<<"c bits:"<<endl;
//    for (int i=0; i<bitsSize; i++){
//        cout<<temp[i];
//    }
//    cout<<endl;

    FieldType two(2);
    FieldType compare = bitwiseLessThan(cBits.data(), rBits.data(), bitsSize);

//    temp[0] = compare;
//    openShare(temp, temp, T);
//     cout<<"c<r? "<<temp[0]<<endl;
    FieldType xorBits0 = cBits[bitsSize - 1] + rBits[bitsSize - 1] - two*cBits[bitsSize - 1]*rBits[bitsSize - 1];
//    temp[0] = xorBits0;
//    openShare(temp, temp, T);
//    cout<<"c0 ^ r0 = "<<temp[0]<<endl;
    vector<FieldType> left(2);
    vector<FieldType> right(2);
    vector<FieldType> mults(2);
    left[0] = compare;
    right[0] = *field->GetOne()-xorBits0;
    left[1] = *field->GetOne() - compare;
    right[1] = xorBits0;


    mults[0] = compare * (*field->GetOne()-xorBits0);
    mults[1] = (*field->GetOne() - compare) * xorBits0;
    DNHonestMultiplication(left.data(), right.data(), mults, 2);
    FieldType res = mults[0] + mults[1];
//    temp[0] = res;
//    openShare(temp, temp, T);
//    cout<<"res = "<<temp[0]<<endl;
    return res;

}

template <class FieldType>
FieldType ProtocolParty<FieldType>::compare(FieldType & a, FieldType & b, int bitSize){


    FieldType two(2);
    //w = [a<p/2]
    vector<FieldType> w(1);
    w[0] = a*two;
    w[0] = *field->GetOne() - LSB(w[0], bitSize);

    //x = [b<p/2]
    vector<FieldType> x(1);
    x[0] = b*two;
    x[0] = *field->GetOne() - LSB(x[0], bitSize);

    //y = [(a-b)<p/2]
    vector<FieldType> y(1);
    y[0] = (a-b)*two;
    y[0] = *field->GetOne() - LSB(y[0], bitSize);

    //compute result = w(x+y-2xy) + 1-y-x+xy

    //xy
    vector<FieldType> xy(1);
    DNHonestMultiplication(x.data(), y.data(), xy, 1);

    //x+y-2xy
    vector<FieldType> calc(1);
    calc[0] = x[0] + y[0] - two*xy[0];

    //w(x+y-2xy)
    DNHonestMultiplication(w.data(), calc.data(), calc, 1);

    //result = w(x+y-2xy) + 1-y-x+xy
    FieldType res = calc[0] + *field->GetOne() - y[0] - x[0] + xy[0];

    return res;

}

template <class FieldType>
void ProtocolParty<FieldType>::getBits(FieldType x, vector<FieldType> & bits, int size) {
    for (int i=size-1; i>=0; i--) {
        bits[i] = x.elem & 1;
        x.elem =  x.elem >> 1;
    }
}

/**
 * some global variables are initialized
 */
template <class FieldType>
void ProtocolParty<FieldType>::initializationPhase() {
//    bigR.resize(1);
//
    beta.resize(1);
    y_for_interpolate.resize(N);

    alpha.resize(N); // N distinct non-zero field elements
//    vector<FieldType> alpha1(N-T);
//    vector<FieldType> alpha2(T);
//
    beta[0] = field->GetElement(0); // zero of the field
    matrix_for_interpolate.allocate(1,N, field);
//
//
    matrix_him.allocate(N,N,field);
    matrix_vand.allocate(N,N,field);
    matrix_vand_transpose.allocate(N,N,field);
//    m.allocate(T, N-T,field);

    // Compute Vandermonde matrix VDM[i,k] = alpha[i]^k
    matrix_vand.InitVDM();
    matrix_vand_transpose.InitVDMTranspose();

    // Prepare an N-by-N hyper-invertible matrix
    matrix_him.InitHIM();

    // N distinct non-zero field elements
    for(int i=0; i<N; i++)
    {
        alpha[i]=field->GetElement(i+1);
    }
//
//    for(int i = 0; i < N-T; i++)
//    {
//        alpha1[i] = alpha[i];
//    }
//    for(int i = N-T; i < N; i++)
//    {
//        alpha2[i - (N-T)] = alpha[i];
//    }
//
//    m.InitHIMByVectors(alpha1, alpha2);
//
    matrix_for_interpolate.InitHIMByVectors(alpha, beta);
//
//    vector<FieldType> alpha_until_t(T + 1);
//    vector<FieldType> alpha_from_t(N - 1 - T);
//
//    // Interpolate first d+1 positions of (alpha,x)
//    matrix_for_t.allocate(N - 1 - T, T + 1, field); // slices, only positions from 0..d
//    //matrix_for_t.InitHIMByVectors(alpha_until_t, alpha_from_t);
//    matrix_for_t.InitHIMVectorAndsizes(alpha, T+1, N-T-1);
//
//    vector<FieldType> alpha_until_2t(2*T + 1);
//    vector<FieldType> alpha_from_2t(N - 1 - 2*T);
//
//    // Interpolate first d+1 positions of (alpha,x)
//    matrix_for_2t.allocate(N - 1 - 2*T, 2*T + 1, field); // slices, only positions from 0..d
//    //matrix_for_2t.InitHIMByVectors(alpha_until_2t, alpha_from_2t);
//    matrix_for_2t.InitHIMVectorAndsizes(alpha, 2*T + 1, N-(2*T +1));
//
//
//    if(flag_print){
//        cout<< "matrix_for_t : " <<endl;
//        matrix_for_t.Print();
//
//        cout<< "matrix_for_2t : " <<endl;
//        matrix_for_2t.Print();
//
//    }




//    auto t1 = high_resolution_clock::now();
//    readclientsinputs(msgsVectorsFlat, squaresVectorsFlat, countersVectorsFlat, unitVectorsFlat);
//
//    auto t2 = high_resolution_clock::now();
//
//    auto duration = duration_cast<milliseconds>(t2-t1).count();
//    if(flag_print_timings) {
//        cout << "time in milliseconds read clients inputs: " << duration << endl;
//    }

    int polySize;
    vector<FieldType> x(delta+1);
    vector<FieldType> y(delta+1);

    polynomials.resize(delta-1);

    for (int i = 0; i < delta+1; i++) {
        x[i] = i + 1;
        y[i] = (i == 0) ? 0 : 1;
    }
    for (int j=2; j<=delta; j++) {
        polySize = j+1;
        polynomials[j-2].resize(polySize);

        cout << "before interpolate" << endl;
        interpolatePolynomial(polynomials[j-2], x.data(), y.data(), polySize);
        cout << "after interpolate" << endl;
    }

    twoSquares.resize(32);
    twoSquares[31] = 1;
    twoSquares[30] = 2;
    for (int i=29; i>=0; i--){
        twoSquares[i] = twoSquares[i+1]*2;
    }

}


template <class FieldType>
bool ProtocolParty<FieldType>::preparationPhase() {

    int keysize = 16/field->getElementSizeInBytes() + 1;

    int numOfRandomShares = 100000*keysize + 1;
    randomSharesArray.resize(numOfRandomShares);

    //generate enough random shares for the AES key
//    generateRandomShares(numOfRandomShares, randomSharesArray);
    RandomSharingForInputs(numOfRandomShares, randomSharesArray);

    //run offline for all the future multiplications including the multiplication of the protocol

    randomTAnd2TSharesOffset = 0;
    generateRandom2TAndTShares(numParties*100000,randomTAnd2TShares);

    randomWithBitsOffset = 0;
    geneateRandomWithBits(field->getElementSizeInBits(), randomsToLSB, bitsToLSB, 3*1000);
//
//    //first generate numOfTriples random shares
//    generateRandomSharesWithCheck(1, bigR);
//
//    //set this random share to an entire array so we can use the semi honest multiplication
//    bigRVec.resize(numClients*securityParamter);
//    fill(bigRVec.begin(), bigRVec.end(), bigR[0]);


    return true;
}


template <class FieldType>
void ProtocolParty<FieldType>::generateRandom2TAndTShares(int numOfRandomPairs, vector<FieldType>& randomElementsToFill){
    auto t1 = high_resolution_clock::now();
    int robin = 0;
    int no_random = numOfRandomPairs;

    // the number of buckets (each bucket requires one double-sharing
    // from each party and gives N-2T random double-sharings)
    int no_buckets = (no_random / (N-T))+1;

    vector<vector<FieldType>> sendBufsElements(N, vector<FieldType>(no_buckets*2));
    vector<vector<FieldType>> recBufsElements(N, vector<FieldType>(no_buckets*2));

    //maybe add some elements if a partial bucket is needed
    randomElementsToFill.resize(no_buckets*(N-T)*2);
    vector<FieldType> randomElementsOnlyTshares (no_buckets*(N-T) );

    int sizeForEachThread;
    if (no_buckets <= numThreads){
        numThreads = no_buckets;
        sizeForEachThread = 1;
    } else{
        sizeForEachThread = (no_buckets + numThreads - 1)/ numThreads;
    }
    vector<thread> threads(numThreads);

    auto t2 = high_resolution_clock::now();

    auto duration = duration_cast<milliseconds>(t2-t1).count();
    if(flag_print_timings) {
        cout << "time in milliseconds generateRandom2TAndTShares: " << duration << endl;
    }
    t1 = high_resolution_clock::now();




    for (int t=0; t<numThreads; t++) {

        if ((t + 1) * sizeForEachThread <= no_buckets) {
            threads[t] = thread(&ProtocolParty::calcSendBufElements, this, ref(sendBufsElements), ref(prgs[t]), t * sizeForEachThread, (t + 1) * sizeForEachThread);
        } else {
            threads[t] = thread(&ProtocolParty::calcSendBufElements, this, ref(sendBufsElements), ref(prgs[t]), t * sizeForEachThread, no_buckets);
        }
    }
    for (int t=0; t<numThreads; t++){
        threads[t].join();
    }

    t2 = high_resolution_clock::now();

    duration = duration_cast<milliseconds>(t2-t1).count();
    if(flag_print_timings) {
        cout << "time in milliseconds calcSendBufElements: " << duration << endl;
    }


    t1 = high_resolution_clock::now();


    roundFunctionSyncElements(sendBufsElements, recBufsElements,4);
    t2 = high_resolution_clock::now();

    duration = duration_cast<milliseconds>(t2-t1).count();
    if(flag_print_timings) {
        cout << "time in milliseconds round function: " << duration << endl;
    }

    t1 = high_resolution_clock::now();

    for (int t=0; t<numThreads; t++) {

        if ((t + 1) * sizeForEachThread <= no_buckets) {
            threads[t] = thread(&ProtocolParty::calcRecBufElements, this, ref(recBufsElements), ref(randomElementsToFill), t * sizeForEachThread, (t + 1) * sizeForEachThread);
        } else {
            threads[t] = thread(&ProtocolParty::calcRecBufElements, this, ref(recBufsElements), ref(randomElementsToFill), t * sizeForEachThread, no_buckets);
        }
    }
    for (int t=0; t<numThreads; t++){
        threads[t].join();
    }

    t2 = high_resolution_clock::now();

    duration = duration_cast<milliseconds>(t2-t1).count();
    if(flag_print_timings) {
        cout << "time in milliseconds calcRecBufElements: " << duration << endl;
    }

    //check validity of the t-shares. 2t-shares do not have to be checked
    //copy the t-shares for checking
    t1 = high_resolution_clock::now();
    for(int i=0; i<randomElementsOnlyTshares.size(); i++){

        randomElementsOnlyTshares[i] = randomElementsToFill[2*i];
    }

    t2 = high_resolution_clock::now();

    duration = duration_cast<milliseconds>(t2-t1).count();
    if(flag_print_timings) {
        cout << "time in milliseconds copy: " << duration << endl;
    }

//    t1 = high_resolution_clock::now();
//    batchConsistencyCheckOfShares(randomElementsOnlyTshares);
//
//    t2 = high_resolution_clock::now();
//
//    duration = duration_cast<milliseconds>(t2-t1).count();
//    if(flag_print_timings) {
//        cout << "time in milliseconds batch consistency: " << duration << endl;
//    }

}


template <class FieldType>
void ProtocolParty<FieldType>::generateRandomShares(int numOfRandoms, vector<FieldType> &randomElementsToFill) {
    int index = 0;
    int no_random = numOfRandoms;

    vector<FieldType> x1(N),y1(N), x2(N),y2(N), t1(N), r1(N), t2(N), r2(N);

    vector<vector<FieldType>> sendBufsElements(N);
    vector<vector<FieldType>> recBufsElements(N);

    // the number of buckets (each bucket requires one double-sharing
    // from each party and gives N-2T random double-sharings)
    int no_buckets = (no_random / (N - T)) + 1;

    //sharingBufTElements.resize(no_buckets*(N-2*T)); // my shares of the double-sharings
    //sharingBuf2TElements.resize(no_buckets*(N-2*T)); // my shares of the double-sharings

    //maybe add some elements if a partial bucket is needed
    randomElementsToFill.resize(no_buckets*(N - T));


    for(int i=0; i < N; i++)
    {
        sendBufsElements[i].resize(no_buckets);
        recBufsElements[i].resize(no_buckets);
    }

    /**
     *  generate random sharings.
     *  first degree t.
     *
     */
    for(int k=0; k < no_buckets; k++)
    {
        // generate random degree-T polynomial
        for(int i = 0; i < T + 1; i++)
        {
            // A random field element, uniform distribution, note that x1[0] is the secret which is also random
            x1[i] = field->Random();

        }

        matrix_vand.MatrixMult(x1, y1, T + 1); // eval poly at alpha-positions

        // prepare shares to be sent
        for(int i=0; i < N; i++)
        {
            //cout << "y1[ " <<i<< "]" <<y1[i] << endl;
            sendBufsElements[i][k] = y1[i];

        }
    }

    if(flag_print) {
        for (int i = 0; i < N; i++) {
            for (int k = 0; k < sendBufsElements[0].size(); k++) {

                // cout << "before roundfunction4 send to " <<i <<" element: "<< k << " " << sendBufsElements[i][k] << endl;
            }
        }
        cout << "sendBufs" << endl;
        cout << "N" << N << endl;
        cout << "T" << T << endl;
    }

    roundFunctionSyncElements(sendBufsElements, recBufsElements, 4);


    for(int k=0; k < no_buckets; k++) {
        for (int i = 0; i < N; i++) {
            t1[i] = recBufsElements[i][k];

        }
        matrix_vand_transpose.MatrixMult(t1, r1, N - T);

        //copy the resulting vector to the array of randoms
        for (int i = 0; i < N - T; i++) {

            randomElementsToFill[index] = r1[i];
            index++;

        }
    }
}


template <class FieldType>
void ProtocolParty<FieldType>::calcSendBufElements(vector<vector<FieldType>> & sendBufsElements, PrgFromOpenSSLAES & prg, int start, int end){

    vector<FieldType> x1(N),y1(N), x2(N),y2(N);

    int* tempInt;

    for(int k=start; k < end; k++)
    {
        // generate random degree-T polynomial
        tempInt = (int*)prg.getPRGBytesEX((T+1)*4);

        for(int i = 0; i < T+1; i++)
        {
            // A random field element, uniform distribution, note that x1[0] is the secret which is also random
            x1[i] = field->GetElement(tempInt[i]);

        }

        matrix_vand.MatrixMult(x1, y1,T+1); // eval poly at alpha-positions

        x2[0] = x1[0];

        // generate random degree-T polynomial
        tempInt = (int*)prg.getPRGBytesEX((2*T+1)*4);

        for(int i = 1; i < 2*T+1; i++)
        {
            // A random field element, uniform distribution, note that x1[0] is the secret which is also random
            x2[i] = field->GetElement(tempInt[i]);

        }

        matrix_vand.MatrixMult(x2, y2,2*T+1);

        // prepare shares to be sent
        for(int i=0; i < N; i++)
        {
            //cout << "y1[ " <<i<< "]" <<y1[i] << endl;
            sendBufsElements[i][2*k] = y1[i];
            sendBufsElements[i][2*k + 1] = y2[i];

        }
    }
}

template <class FieldType>
void ProtocolParty<FieldType>::calcRecBufElements(vector<vector<FieldType>> & recBufsElements, vector<FieldType> & randomElementsToFill, int start, int end){

    vector<FieldType> t1(N), r1(N), t2(N), r2(N);

    for(int k=start; k < end; k++) {
        for (int i = 0; i < N; i++) {

            t1[i] = recBufsElements[i][2*k];
            t2[i] = recBufsElements[i][(2*k +1)];


        }
        matrix_vand_transpose.MatrixMult(t1, r1,N-T);
        matrix_vand_transpose.MatrixMult(t2, r2,N-T);

        //copy the resulting vector to the array of randoms
        for (int i = 0; i < (N - T); i++) {

            randomElementsToFill[k*2] = r1[i];
            randomElementsToFill[k*2 +1] = r2[i];

        }
    }
}


template <class FieldType>
void ProtocolParty<FieldType>::DNHonestMultiplication(FieldType *a, FieldType *b, vector<FieldType> &cToFill, int numOfTupples) {

    int fieldByteSize = field->getElementSizeInBytes();
    vector<FieldType> xyMinusRShares(numOfTupples);//hold both in the same vector to send in one batch
    vector<byte> xyMinusRSharesBytes(numOfTupples *fieldByteSize);//hold both in the same vector to send in one batch

    vector<FieldType> xyMinusR;//hold both in the same vector to send in one batch
    vector<byte> xyMinusRBytes;

    vector<vector<FieldType>> recBufsElements(N);
    vector<vector<FieldType>> sendBufsElements(N);


    //generate the shares for x+a and y+b. do it in the same array to send once
    for (int k = 0; k < numOfTupples; k++)//go over only the logit gates
    {
        //compute the share of xy-r
        xyMinusRShares[k] = a[k]*b[k] - randomTAnd2TShares[randomTAnd2TSharesOffset + 2*k+1];

    }

    //set the acctual number of mult gate proccessed in this layer
    int acctualNumOfMultGates = numOfTupples;
    int numOfElementsForParties = acctualNumOfMultGates/N;
    int indexForDecreasingSize = acctualNumOfMultGates - numOfElementsForParties *N;

    int counter=0;
    int currentNumOfElements;
    for(int i=0; i<N; i++){

        currentNumOfElements = numOfElementsForParties;
        if(i<indexForDecreasingSize)
            currentNumOfElements++;

        //fill the send buf according to the number of elements to send to each party
        sendBufsElements[i].resize(currentNumOfElements);

        for(int j=0; j<currentNumOfElements; j++) {

            sendBufsElements[i][j] = xyMinusRShares[counter];
            counter++;

        }

    }

    //resize the recbuf array.
    int myNumOfElementsToExpect = numOfElementsForParties;
    if (partyId < indexForDecreasingSize) {
        myNumOfElementsToExpect = numOfElementsForParties + 1;
    }

    for(int i=0;i<N;i++){

        //recBufsBytes[i].resize(myNumOfElementsToExpect*fieldByteSize);
        recBufsElements[i].resize(myNumOfElementsToExpect);

    }


    roundFunctionSyncElements(sendBufsElements, recBufsElements,20);

    xyMinusR.resize(myNumOfElementsToExpect);
    xyMinusRBytes.resize(myNumOfElementsToExpect*fieldByteSize);

    //reconstruct the shares that I am responsible of recieved from the other parties
    vector<FieldType> xyMinurAllShares(N);

    for (int k = 0;k < myNumOfElementsToExpect; k++)//go over only the logit gates
    {
        for (int i = 0; i < N; i++) {

            xyMinurAllShares[i] = recBufsElements[i][k];
        }

        // reconstruct the shares by P0
        xyMinusR[k] = interpolate(xyMinurAllShares);

    }


    //prepare the send buffers
    for(int i=0; i<N; i++){
        //sendBufsBytes[i] = xyMinusRBytes;
        sendBufsElements[i] = xyMinusR;
    }


    for(int i=0; i<N; i++){

        currentNumOfElements = numOfElementsForParties;
        if(i<indexForDecreasingSize)
            currentNumOfElements++;

        //recBufsBytes[i].resize(currentNumOfElements* fieldByteSize);
        recBufsElements[i].resize(currentNumOfElements);

    }
    roundFunctionSyncElements(sendBufsElements, recBufsElements,21);


    xyMinusR.resize(acctualNumOfMultGates);
    counter = 0;

    for(int i=0; i<N; i++){

        currentNumOfElements = numOfElementsForParties;
        if(i<indexForDecreasingSize)
            currentNumOfElements++;

        //fill the send buf according to the number of elements to send to each party
        for(int j=0; j<currentNumOfElements; j++) {

            //xyMinusR[counter] = field->bytesToElement(recBufsBytes[i].data() + (j * fieldByteSize));
            xyMinusR[counter] = recBufsElements[i][j];

            counter++;

        }

    }


    for (int k = 0; k < numOfTupples; k++)
    {
        cToFill[k] = randomTAnd2TShares[randomTAnd2TSharesOffset + 2*k] + xyMinusR[k];
    }

    randomTAnd2TSharesOffset+=numOfTupples*2;

}

template <class FieldType>
void ProtocolParty<FieldType>::openShare(vector<FieldType> &shares, vector<FieldType> &secrets, int d){


    vector<vector<byte>> recBufsBytes(N);

    openShareSetRecBuf(shares, secrets, d, recBufsBytes);

}
template <class FieldType>
void ProtocolParty<FieldType>::openShareSetRecBuf(vector<FieldType> &shares, vector<FieldType> &secrets,
                                                  int d, vector<vector<byte>> &recBufsBytes){

    int numOfRandomShares = shares.size();
    vector<vector<byte>> sendBufsBytes(N);

    vector<FieldType> x1(N);
    int fieldByteSize = field->getElementSizeInBytes();

    //calc the number of elements needed for 128 bit AES key

    //resize vectors
    for(int i=0; i < N; i++)
    {
        sendBufsBytes[i].resize(numOfRandomShares*fieldByteSize);
        recBufsBytes[i].resize(numOfRandomShares*fieldByteSize);
    }

    //set the first sending data buffer
    for(int j=0; j<numOfRandomShares;j++) {
        field->elementToBytes(sendBufsBytes[0].data() + (j * fieldByteSize), shares[j]);
    }

    //copy the same data for all parties
    for(int i=1; i<N; i++){

        sendBufsBytes[i] = sendBufsBytes[0];
    }

    //call the round function to send the shares to all the users and get the other parties share
    roundFunctionSync(sendBufsBytes, recBufsBytes,12);

    //reconstruct each set of shares to get the secret

    for(int k=0; k<numOfRandomShares; k++){

        //get the set of shares for each element
        for(int i=0; i < N; i++) {

            x1[i] = field->bytesToElement(recBufsBytes[i].data() + (k*fieldByteSize));
        }


        secrets[k] = reconstructShare(x1, d);

    }

}


template <class FieldType>
void ProtocolParty<FieldType>::roundFunctionSync(vector<vector<byte>> &sendBufs, vector<vector<byte>> &recBufs, int round) {

    //cout<<"in roundFunctionSync "<< round<< endl;

    int numThreads = parties.size();
    int numPartiesForEachThread;

    if (parties.size() <= numThreads){
        numThreads = parties.size();
        numPartiesForEachThread = 1;
    } else{
        numPartiesForEachThread = (parties.size() + numThreads - 1)/ numThreads;
    }


    recBufs[partyId] = move(sendBufs[partyId]);
    //recieve the data using threads
    vector<thread> threads(numThreads);
    for (int t=0; t<numThreads; t++) {
        if ((t + 1) * numPartiesForEachThread <= parties.size()) {
            threads[t] = thread(&ProtocolParty::exchangeData, this, ref(sendBufs), ref(recBufs),
                                t * numPartiesForEachThread, (t + 1) * numPartiesForEachThread);
        } else {
            threads[t] = thread(&ProtocolParty::exchangeData, this, ref(sendBufs), ref(recBufs), t * numPartiesForEachThread, parties.size());
        }
    }
    for (int t=0; t<numThreads; t++){
        threads[t].join();
    }

}


template <class FieldType>
void ProtocolParty<FieldType>::exchangeData(vector<vector<byte>> &sendBufs, vector<vector<byte>> &recBufs, int first, int last){


    //cout<<"in exchangeData";
    for (int i=first; i < last; i++) {

        if ((partyId) < parties[i]->getID()) {


            //send shares to my input bits
            parties[i]->getChannel()->write(sendBufs[parties[i]->getID()].data(), sendBufs[parties[i]->getID()].size());
            //cout<<"write the data:: my Id = " << m_partyId - 1<< "other ID = "<< parties[i]->getID() <<endl;


            //receive shares from the other party and set them in the shares array
            parties[i]->getChannel()->read(recBufs[parties[i]->getID()].data(), recBufs[parties[i]->getID()].size());
            //cout<<"read the data:: my Id = " << m_partyId-1<< "other ID = "<< parties[i]->getID()<<endl;

        } else{


            //receive shares from the other party and set them in the shares array
            parties[i]->getChannel()->read(recBufs[parties[i]->getID()].data(), recBufs[parties[i]->getID()].size());
            //cout<<"read the data:: my Id = " << m_partyId-1<< "other ID = "<< parties[i]->getID()<<endl;



            //send shares to my input bits
            parties[i]->getChannel()->write(sendBufs[parties[i]->getID()].data(), sendBufs[parties[i]->getID()].size());
            //cout<<"write the data:: my Id = " << m_partyId-1<< "other ID = "<< parties[i]->getID() <<endl;


        }

    }


}


template <class FieldType>
void ProtocolParty<FieldType>::roundFunctionSyncElements(vector<vector<FieldType>> &sendBufs, vector<vector<FieldType>> &recBufs, int round) {

    //cout<<"in roundFunctionSync "<< round<< endl;

    int numThreads = parties.size();
    int numPartiesForEachThread;

    if (parties.size() <= numThreads){
        numThreads = parties.size();
        numPartiesForEachThread = 1;
    } else{
        numPartiesForEachThread = (parties.size() + numThreads - 1)/ numThreads;
    }


    recBufs[partyId] = move(sendBufs[partyId]);
    //recieve the data using threads
    vector<thread> threads(numThreads);
    for (int t=0; t<numThreads; t++) {
        if ((t + 1) * numPartiesForEachThread <= parties.size()) {
            threads[t] = thread(&ProtocolParty::exchangeDataElements, this, ref(sendBufs), ref(recBufs),
                                t * numPartiesForEachThread, (t + 1) * numPartiesForEachThread);
        } else {
            threads[t] = thread(&ProtocolParty::exchangeDataElements, this, ref(sendBufs), ref(recBufs), t * numPartiesForEachThread, parties.size());
        }
    }
    for (int t=0; t<numThreads; t++){
        threads[t].join();
    }

}


template <class FieldType>
void ProtocolParty<FieldType>::exchangeDataElements(vector<vector<FieldType>> &sendBufs, vector<vector<FieldType>> &recBufs, int first, int last) {


    //cout<<"in exchangeData";
    for (int i = first; i < last; i++) {

        if ((partyId) < parties[i]->getID()) {


            if (sendBufs[parties[i]->getID()].size() > 0) {
                //send shares to my input bits
                parties[i]->getChannel()->write((byte *) sendBufs[parties[i]->getID()].data(),
                                                sendBufs[parties[i]->getID()].size() * field->getElementSizeInBytes());
                //cout<<"write the data:: my Id = " << m_partyId - 1<< "other ID = "<< parties[i]->getID() <<endl;
            }

            if (recBufs[parties[i]->getID()].size() > 0) {
                //receive shares from the other party and set them in the shares array
                parties[i]->getChannel()->read((byte *) recBufs[parties[i]->getID()].data(),
                                               recBufs[parties[i]->getID()].size() * field->getElementSizeInBytes());
                //cout<<"read the data:: my Id = " << m_partyId-1<< "other ID = "<< parties[i]->getID()<<endl;
            }

        } else {

            if (recBufs[parties[i]->getID()].size() > 0) {
                //receive shares from the other party and set them in the shares array
                parties[i]->getChannel()->read((byte *) recBufs[parties[i]->getID()].data(),
                                               recBufs[parties[i]->getID()].size() * field->getElementSizeInBytes());
                //cout<<"read the data:: my Id = " << m_partyId-1<< "other ID = "<< parties[i]->getID()<<endl;
            }

            if (sendBufs[parties[i]->getID()].size() > 0) {

                //send shares to my input bits
                parties[i]->getChannel()->write((byte *) sendBufs[parties[i]->getID()].data(),
                                                sendBufs[parties[i]->getID()].size() * field->getElementSizeInBytes());
                //cout<<"write the data:: my Id = " << m_partyId-1<< "other ID = "<< parties[i]->getID() <<endl;
            }

        }

    }
}


template <class FieldType>
FieldType ProtocolParty<FieldType>::reconstructShare(vector<FieldType>& x, int d){

//    if (!checkConsistency(x, d))
//    {
//        // someone cheated!
//
//        cout << "cheating reconstruct!!!" << '\n';
//        //exit(0);
//    }
//    else
        return interpolate(x);
}


// Interpolate polynomial at position Zero
template <class FieldType>
FieldType ProtocolParty<FieldType>::interpolate(vector<FieldType>& x)
{
    matrix_for_interpolate.MatrixMult(x, y_for_interpolate);
    return y_for_interpolate[0];
}

template <class FieldType>
void ProtocolParty<FieldType>::interpolatePolynomial(vector<FieldType>& coeff, const FieldType* a, FieldType* b, int size)
{
    int m = size;

    FieldType one(1);
    FieldType zero;

    FieldType p(FieldType::p);

    vector<FieldType> prod(size);
    memcpy(prod.data(), a, size*field->getElementSizeInBytes());

    FieldType t1, t2;

    long k, i;

    vector<FieldType> res;
    res.resize(m);

    for (k = 0; k < m; k++) {

        const FieldType& aa = a[k];

        t1 = 1;
        for (i = k-1; i >= 0; i--) {
            t1 = t1*aa; //mul(t1, t1, aa);
            t1 = t1 + prod[i];//add(t1, t1, prod[i]);
        }

        t2 = 0; //clear(t2);
        for (i = k-1; i >= 0; i--) {
            t2 = t2*aa; //mul(t2, t2, aa);
            t2 = t2 + res[i];//add(t2, t2, res[i]);
        }


        t1 = one/t1;//inv(t1, t1);
        t2 = b[k] - t2;//sub(t2, b[k], t2);
        t1 = t1*t2;//mul(t1, t1, t2);

        for (i = 0; i < k; i++) {
            t2 = prod[i]*t1;//mul(t2, prod[i], t1);
            res[i] = res[i] + t2;//add(res[i], res[i], t2);
        }

        res[k] = t1;

        if (k < m-1) {
            if (k == 0)
                prod[0] = p - prod[0];//sub(prod[0], to_ZZ_p(ZZ_pInfo->p),prod[0]);//sub(prod[0], ZZ_p::modulus(), prod[0]);//negate(prod[0], prod[0]);
            else {
                t1 = p - a[k];//sub(t1, to_ZZ_p(ZZ_pInfo->p),a[k]);//negate(t1, a[k]);
                prod[k] = t1 + prod[k-1];//add(prod[k], t1, prod[k-1]);
                for (i = k-1; i >= 1; i--) {
                    t2 = prod[i]*t1;//mul(t2, prod[i], t1);
                    prod[i] = t2 + prod[i-1];//add(prod[i], t2, prod[i-1]);
                }
                prod[0] = prod[0]*t1;//mul(prod[0], prod[0], t1);
            }
        }
    }

    while (m > 0 && !(res[m-1]!=zero)) m--;
    res.resize(m);


    coeff = res;
}

template <class FieldType>
FieldType ProtocolParty<FieldType>::evalPolynomial(const vector<FieldType>& coeff, vector<FieldType> & a)
// does a Horner evaluation
{

//    cout<<"poly:"<<endl;
//    for (int i=0; i<coeff.size(); i++){
//        cout<<coeff[i]<<" ";
//    }
//    cout<<endl;


    int size = coeff.size();
    FieldType acc = coeff[0];
    for (int j=0; j<size - 1; j++) {
        acc += a[j] * coeff[j+1];
//        acc = coeff[coeff.size() - 1];
//
//        for (int i = coeff.size() - 2; i >= 0; i--) {
//            acc = acc * a[j];//mul(acc, acc, a);
//            acc = acc + coeff[i];//add(acc, acc, f.rep[i]);
//        }
//
//        b[j] = acc;
    }

    return acc;
}


template <class FieldType>
bool ProtocolParty<FieldType>::RandomSharingForInputs(int no_random, vector<FieldType> & sharingBufInputsTElements)
{
    vector<vector<byte>> recBufsBytes(N);
    //vector<vector<byte>> recBufs1Bytes(N);
    int robin = 0;

    // the number of random double sharings we need altogether
    vector<FieldType> x1(N),y1(N);

    vector<vector<FieldType>> sendBufsElements(N);
    vector<vector<byte>> sendBufsBytes(N);

    // the number of buckets (each bucket requires one double-sharing
    // from each party and gives N-2T random double-sharings)
    int no_buckets = (no_random / (N-2*T))+1;

    //sharingBufTElements.resize(no_buckets*(N-2*T)); // my shares of the double-sharings
    //sharingBuf2TElements.resize(no_buckets*(N-2*T)); // my shares of the double-sharings
    sharingBufInputsTElements.resize(no_buckets*(N-2*T));

    for(int i=0; i < N; i++)
    {
        //sendBufs[i] = "";

        sendBufsElements[i].resize(no_buckets);
        sendBufsBytes[i].resize(no_buckets*field->getElementSizeInBytes());
        recBufsBytes[i].resize(no_buckets*field->getElementSizeInBytes());
    }

    /**
     *  generate double sharings.
     *  first degree t.
     *  subsequent: degree 2t with same secret.
     */
    high_resolution_clock::time_point t1 = high_resolution_clock::now();
    for(int k=0; k < no_buckets; k++)
    {
        // generate random degree-T polynomial
        for(int i = 0; i < T+1; i++)
        {
            // A random field element, uniform distribution
            x1[i] = field->Random();

        }

        matrix_vand.MatrixMult(x1, y1, T+1); // eval poly at alpha-positions


        // prepare shares to be sent
        for(int i=0; i < N; i++)
        {
            //cout << "y1[ " <<i<< "]" <<y1[i] << endl;
            sendBufsElements[i][k] = y1[i];

        }




    }//end print one



    if(flag_print) {
        for (int i = 0; i < N; i++) {
            for (int k = 0; k < sendBufsElements[0].size(); k++) {

                // cout << "before roundfunction4 send to " <<i <<" element: "<< k << " " << sendBufsElements[i][k] << endl;
            }
        }
    }

    high_resolution_clock::time_point t2 = high_resolution_clock::now();
    auto duration = duration_cast<milliseconds>( t2 - t1 ).count();
    //cout << "generate random degree-T polynomial took : " <<duration<<" ms"<<endl;

    if(flag_print) {
        cout << "sendBufs" << endl;
        cout << "N" << N << endl;
        cout << "T" << T << endl;



    }
//
    int fieldByteSize = field->getElementSizeInBytes();
    for(int i=0; i < N; i++)
    {
        for(int j=0; j<sendBufsElements[i].size();j++) {
            field->elementToBytes(sendBufsBytes[i].data() + (j * fieldByteSize), sendBufsElements[i][j]);
        }
    }


    high_resolution_clock::time_point t3 = high_resolution_clock::now();
    //comm->roundfunctionI(sendBufsBytes, recBufsBytes,4);
    roundFunctionSync(sendBufsBytes, recBufsBytes,4);
    high_resolution_clock::time_point t4 = high_resolution_clock::now();
    auto duration2 = duration_cast<milliseconds>( t4 - t3 ).count();
    //cout << "roundfunctionI took : " <<duration2<<" ms"<<endl;






    if(flag_print) {
        for (int i = 0; i < N; i++) {
            for (int k = 0; k < sendBufsBytes[0].size(); k++) {

                cout << "roundfunction4 send to " <<i <<" element: "<< k << " " << (int)sendBufsBytes[i][k] << endl;
            }
        }
    }

    //cout << endl;

    if(flag_print) {
        for (int i = 0; i < N; i++) {
            for (int k = 0; k < recBufsBytes[0].size(); k++) {
                cout << "roundfunction4 receive from " <<i <<" element: "<< k << " " << (int) recBufsBytes[i][k] << endl;
            }
        }
    }


    /**
     * Apply hyper-invertible matrix on each bucket.
     * From the resulting sharings, 2T are being reconstructed towards some party,
     * the remaining N-2T are kept as prepared sharings.
     * For balancing, we do round-robin the party how shall reconstruct and check!
     */


    //vector<vector<FieldType>> sendBufs1Elements(N);
    //vector<vector<byte>> sendBufs1Bytes(N);

    for(int i=0; i<N; i++){
        sendBufsElements[i].clear();

    }

    int fieldBytesSize = field->getElementSizeInBytes();

    // x1 : used for the N degree-t sharings
    // x2 : used for the N degree-2t sharings
    for(int k=0; k < no_buckets; k++) {
        // generate random degree-T polynomial
        for (int i = 0; i < N; i++) {
            x1[i] = field->bytesToElement(recBufsBytes[i].data() + (k*fieldBytesSize));

        }
        matrix_him.MatrixMult(x1, y1);
        // these shall be checked
        for (int i = 0; i < 2 * T; i++) {
            sendBufsElements[robin].push_back(y1[i]);
            robin = (robin+1) % N; // next robin

        }
        // Y1 : the degree-t shares of my poly
        // Y2 : the degree 2t shares of my poly
        for (int i = 2 * T; i < N; i++) {

            sharingBufInputsTElements[k*(N-2*T) + i - 2*T] = y1[i];
            //sharingBufTElements[k*(N-2*T) + i - 2*T] = y1[i];
            //sharingBuf2TElements[k*(N-2*T) + i - 2*T] =  y2[i];
        }

    }

    for(int i=0; i < N; i++)
    {
        sendBufsBytes[i].resize(sendBufsElements[i].size()*fieldByteSize);
        //cout<< "size of sendBufs1Elements["<<i<<" ].size() is " << sendBufs1Elements[i].size() <<"myID =" <<  m_partyId<<endl;
        recBufsBytes[i].resize(sendBufsElements[partyId].size()*fieldByteSize);
        for(int j=0; j<sendBufsElements[i].size();j++) {
            field->elementToBytes(sendBufsBytes[i].data() + (j * fieldByteSize), sendBufsElements[i][j]);
        }
    }


    if(flag_print)
        cout << "before round" << endl;

    t3 = high_resolution_clock::now();
    roundFunctionSync(sendBufsBytes, recBufsBytes,5);
    t4 = high_resolution_clock::now();
    duration2 = duration_cast<milliseconds>( t4 - t3 ).count();

    if(flag_print) {
        cout << "after round" << endl;}
    int count = no_buckets * (2*T) / N; // nr of sharings *I* have to check
    // got one in the last round
    if(no_buckets * (2*T)%N > partyId) { // maybe -1
        count++;
    }


    for(int k=0; k < count; k++) {
        for (int i = 0; i < N; i++) {

            x1[i] = field->bytesToElement(recBufsBytes[i].data() + (k*fieldBytesSize));
        }


        vector<FieldType> x_until_d(N);
        for(int i=0; i<T; i++)
        {
            x_until_d[i] = x1[i];
        }
        for(int i=T; i<N; i++)
        {
            x_until_d[i] = *(field->GetZero());
        }
        if(flag_print) {
            cout << "k " << k << "interpolate(x1).toString()  " << field->elementToString(interpolate(x1)) << endl;
        }
    }
    return true;
}


#endif //COMPARISON_PROTOCOLPARTY_H
