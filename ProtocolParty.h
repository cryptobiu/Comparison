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
    int randomTAnd2TSharesOffset;
    vector<FieldType> randomSharesArray;
    int randomSharesOffset;

    vector<FieldType> randomsToLSB;
    vector<vector<FieldType>> bitsToLSB;
    int randomWithBitsOffset;

    int numCompares;
    int debugDuration;


    void initializationPhase();
    bool preparationPhase();
    void interpolatePolynomial(vector<FieldType>& coeff, const FieldType* a, FieldType* b, int size);
    FieldType evalPolynomial(const vector<FieldType>& coeff, FieldType* a);

    bool RandomSharingForInputs(int no_random, vector<FieldType> & sharingBufInputsTElements);
    void generateRandom2TAndTShares(int numOfRandomPairs, vector<FieldType>& randomElementsToFill);
    void generateRandomShares(int numOfRandoms, vector<FieldType> &randomElementsToFill);

    void calcSendBufElements(vector<vector<FieldType>> & sendBufsElements, PrgFromOpenSSLAES & prg, int start, int end);
    void calcRecBufElements(vector<vector<FieldType>> & recBufsElements, vector<FieldType> & randomElementsToFill, int start, int end);

    void openShare(vector<FieldType> &shares, vector<FieldType> &secrets, int d);

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

    void DNHonestMultiplication(FieldType *a, FieldType *b, vector<FieldType> &cToFill, int numOfTupples, bool sop = false);
    void inverse(FieldType *a, vector<FieldType> &bToFill);

    void randomBit(vector<FieldType> &bToFill);
    void checkRandomBit(vector<FieldType> &bits);

    void unboundedOR(vector<vector<FieldType>> & shares, int sizeTotal, vector<FieldType> & resToFill, int offset);
    void checkUnboundedOR(FieldType* bits, int numBits);

//    void prefixOR(FieldType* shares, vector<FieldType> & bToFill);
    void prefixOR(vector<vector<FieldType>> & shares, vector<vector<FieldType>> & bToFill);

//    FieldType bitwiseLessThan(FieldType* a, FieldType* b, int size);
    void bitwiseLessThan(vector<vector<FieldType>> & a, vector<vector<FieldType>> & b, vector<FieldType> & result);

    void geneateRandomWithBits(int size, vector<FieldType> & randomToFill, vector<vector<FieldType>> & bitsToFill, int numRandoms);
    void getBits(vector<FieldType> & x, vector<vector<FieldType>> & bits, int size);
//    FieldType LSB(FieldType & a, int bitsSize);
    void LSB(vector<FieldType> & x, int bitsSize, vector<FieldType> & lsb);

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
    } else if(fieldType.compare("ZpMersenne13") == 0) {
        field = new TemplateField<FieldType>(8191);
    }
    else if(fieldType.compare("Zp16BitPrime") == 0) {
        field = new TemplateField<FieldType>(8191);
    }

    N = numParties;
    T = (numParties+1)/2 - 1;

    delta = stoi(this->getParser().getValueByKey(arguments, "delta"));
    lambda = stoi(this->getParser().getValueByKey(arguments, "lambda"));
//    delta = 5;
//    lambda = 3;

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

    string tmp = "init times";
    //cout<<"before sending any data"<<endl;
    byte tmpBytes[20];
    for (int i=0; i<parties.size(); i++){
        if (parties[i]->getID() < partyId){
            parties[i]->getChannel()->write(tmp);
            parties[i]->getChannel()->read(tmpBytes, tmp.size());
        } else {
            parties[i]->getChannel()->read(tmpBytes, tmp.size());
            parties[i]->getChannel()->write(tmp);
        }
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
        cout << "time in milliseconds comparePhase: " << duration << endl;
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
    auto start = high_resolution_clock::now();
    for (int i=0; i<numCompares; i++) {

        aC = randomSharesArray[randomSharesOffset++];
        bC = randomSharesArray[randomSharesOffset++];

        res = compare(aC, bC, field->getElementSizeInBits());

    }
    auto end = high_resolution_clock::now();
    auto duration = duration_cast<milliseconds>(end- start).count();
//        cout << "compare took " << duration << " ms"<<endl;

    vector<FieldType> temp(3);
    temp[0] = aC;
    temp[1] = bC;
    temp[2] = res;
    openShare(temp, temp, T);
    cout<<"a = "<<temp[0]<<endl;
    cout<<"b = "<<temp[1]<<endl;
    cout<<"a < b ? = "<<temp[2]<<endl;
    cout<<"all prefix or took = "<<debugDuration<<endl;
    cout << "compute " <<numCompares<<" compares took " << duration/numCompares << " ms in average"<<endl;

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
void ProtocolParty<FieldType>::unboundedOR(vector<vector<FieldType>> & shares, int sizeTotal, vector<FieldType> & resToFill, int offset) {
    int numOfGroups = shares.size();
//    cout<<"num of groups in unbounded = "<<numOfGroups<<endl;
    vector<FieldType> A(sizeTotal);
    vector<FieldType> b(sizeTotal);

    vector<int> sizes(numOfGroups);
    vector<int> offsets(numOfGroups);
    int index = 0;
    for (int j=0; j<numOfGroups; j++) {
        sizes[j] = shares[j].size();
        offsets[j] = index;

        A[offsets[j]] = *field->GetOne();
        for (int i = 0; i < sizes[j]; i++) {
            A[offsets[j]] += shares[j][i];
        }
        index += sizes[j];
    }

//    cout<<"sizes:"<<endl;
//    for (int j=0; j<numOfGroups; j++) {
//        cout<<sizes[j]<<" ";
//    }
//    cout<<endl;
//    cout<<"offsets:"<<endl;
//    for (int j=0; j<numOfGroups; j++) {
//        cout<<offsets[j]<<" ";
//    }
//    cout<<endl;

    int originalOffset = randomSharesOffset;
    vector<FieldType> lastB(numOfGroups);
    for (int j=0; j<numOfGroups; j++) {
        b[offsets[j]] = *field->GetOne();
        for (int i = 1; i < sizes[j]; i++) {
            A[offsets[j] + i] = A[offsets[j]];
            b[offsets[j] + i] = randomSharesArray[randomSharesOffset++];
        }
        lastB[j] = randomSharesArray[randomSharesOffset++];
    }

    //interpolate done in initialization phase
    vector<FieldType> randomsInverse(sizeTotal);
    inverse(randomSharesArray.data() + originalOffset, randomsInverse);
//cout<<"after inverse"<<endl;
    vector<FieldType> c(sizeTotal);
    //calc first column - A * all B values
    DNHonestMultiplication(A.data(), b.data(), c, sizeTotal);

    for (int j=0; j<numOfGroups; j++) {
        c[offsets[j]] = A[offsets[j]]; //the first value should be fixed
    }

    //calc second column - A * all B values
//    DNHonestMultiplication(c.data(), randomsInverse.data(), c, l);
    for (int j=0; j<numOfGroups; j++) {
        for (int i = 0; i < sizes[j]; i++) {
            c[offsets[j] + i] = c[offsets[j] + i] * randomsInverse[offsets[j] + i];
        }
    }

    //open C
    vector<FieldType> cOpened(sizeTotal);
    openShare(c, cOpened, 2*T);

//    cout<<"opened c:"<<endl;
//    for (int i=0; i<sizeTotal; i++){
//        cout<<cOpened[i]<<" ";
//    }
//    cout<<endl;

    FieldType cMult;
    //compute all A squares
    for (int j=0; j<numOfGroups; j++) {
        cMult = cOpened[offsets[j]];
        for (int i = 1; i < sizes[j] - 1; i++) {
            cMult *= cOpened[offsets[j] + i];
            A[offsets[j] + i] = cMult * b[offsets[j] + i + 1];
        }
        A[offsets[j] + sizes[j] - 1] = cMult * cOpened[offsets[j] + sizes[j] - 1] * lastB[j];
    }

//    openShare(A, cOpened, T);
//    cout<<"opened A:"<<endl;
//    for (int i=0; i<l; i++){
//        cout<<cOpened[i]<<" ";
//    }
//    cout<<endl;

    for (int j=0; j<numOfGroups; j++) {
        if (sizes[j] == 1){
            resToFill[offset + j] = shares[j][0];
        }
         else resToFill[offset + j] = evalPolynomial(polynomials[sizes[j]-2], A.data() + offsets[j]);
    }
}


template <class FieldType>
void ProtocolParty<FieldType>::prefixOR(vector<vector<FieldType>> & shares, vector<vector<FieldType>> & bToFill) {

    int numVectors = shares.size();
//    cout<<"original size = "<<size<<endl;

    //calculate xi = unbounded or of each delta bits
    vector<FieldType> x(numVectors*lambda);
    vector<FieldType> y(numVectors*lambda);
    vector<FieldType> f(numVectors*lambda);
    vector<FieldType> s(numVectors*lambda);

    vector<vector<FieldType>> sharesVector(numVectors*lambda);

    int totalSize = 0;
    for (int j=0; j<numVectors; j++) {

        //copy the first lambda-1 groups
        for (int i = 0; i < lambda - 1; i++) {
            sharesVector[j*lambda + i].resize(delta);
            memcpy((byte *) sharesVector[j*lambda + i].data(), (byte *) &shares[j][i * delta],
                   delta * field->getElementSizeInBytes());
        }
        //copy the last gropu. Notice that it can contain less than delta elements
        int remain = shares[j].size() % delta;
        if (remain == 0) {
            sharesVector[j*lambda + lambda - 1].resize(delta);
            memcpy((byte *) sharesVector[j*lambda + lambda - 1].data(), (byte *) &shares[j][(lambda - 1) * delta],
                   delta * field->getElementSizeInBytes());
            totalSize += lambda * delta;

        } else {
            sharesVector[j*lambda + lambda - 1].resize(remain);
            memcpy((byte *) sharesVector[j*lambda + lambda - 1].data(), (byte *) &shares[j][(lambda - 1) * delta],
                   remain * field->getElementSizeInBytes());
            totalSize += (lambda - 1) * delta + remain;
        }
    }

//    vector<FieldType> temp;
//    cout<<"sharesAlligned:"<<endl;
//    for (int i=0; i<lambda; i++){
//        temp.resize(sharesVector[i].size());
//        openShare(sharesVector[i],temp,T);
//        for (int i=0; i<temp.size(); i++){
//            cout<<temp[i]<<" ";
//        }
//    }
//    cout<<endl;
//cout<<"before unbounded"<<endl;
//    cout<<"totalSize = "<<totalSize<<endl;
    unboundedOR(sharesVector, totalSize, x, 0);
//    for (int i=0; i<lambda; i++){
//        x[i] = unboundedOR(sharesAlligned + i*delta, delta);
//    }

//    vector<FieldType> temp;
//    temp.resize(x.size());
//    openShare(x,temp,T);
//    cout<<"x:"<<endl;
//    for (int i=0; i<temp.size(); i++){
//        cout<<temp[i]<<" ";
//    }
//    cout<<endl;

    //calculate yi = prefix or of x
    totalSize = 0;
    sharesVector.resize(numVectors*(lambda-1));
    vector<FieldType> tempToOR(numVectors*(lambda-1));
    for (int j=0; j<numVectors; j++) {
        for (int i = 0; i < lambda - 1; i++) {
            sharesVector[j*(lambda-1) + i].resize(i + 2);
            memcpy((byte *) sharesVector[j*(lambda-1) + i].data(), (byte *) &x[j*lambda], (i + 2) * field->getElementSizeInBytes());
            totalSize += i + 2;
        }
        y[j*lambda] = x[j*lambda];
    }

    unboundedOR(sharesVector, totalSize, tempToOR, 0);

    for (int j=0; j<numVectors; j++) {
        for (int i = 0; i < lambda - 1; i++) {
            y[j*lambda + i + 1] = tempToOR[j*(lambda-1) + i];
        }
    }

//    for (int i=1; i<lambda; i++){
//        y[i] = unboundedOR(x.data(), i+1);
//    }

//    openShare(y,temp,T);
//    cout<<"y:"<<endl;
//    for (int i=0; i<temp.size(); i++){
//        cout<<temp[i]<<" ";
//    }
//    cout<<endl;

    //calculate fi = find the group with the first "1" bit
    for (int j=0; j<numVectors; j++) {
        f[j*lambda] = y[j*lambda];
        for (int i = 1; i < lambda; i++) {
            f[j*lambda + i] = y[j*lambda + i] - y[j*lambda + i - 1];
        }
    }

//    openShare(f,temp,T);
//    cout<<"f:"<<endl;
//    for (int i=0; i<temp.size(); i++){
//        cout<<temp[i]<<" ";
//    }
//    cout<<endl;

    //find the real b values for the group which contains the first "1" bit.
    //Note that the group index is unknown
    vector<FieldType> a(numVectors * delta, *field->GetZero());

    for (int k=0; k<numVectors;k++) {
        for (int j = 0; j < delta; j++) {

            for (int i = 0; i < lambda - 1; i++) {
                a[k*delta + j] += f[k*lambda + i] * shares[k][i * delta + j];

            }

            int remain = shares[k].size() % delta;
            if (j < remain) {
                a[k*delta + j] += f[k*lambda + lambda - 1] * shares[k][(lambda - 1) * delta + j];
            }
        }
    }

    DNHonestMultiplication(a.data(), a.data(), a, numVectors * delta, true);

    //calculate the real b values for the group which contains the first "1" bit.
    //Note that the group index is unknown
    vector<FieldType> b(numVectors * delta);

    totalSize = 0;
    sharesVector.resize(numVectors * (delta-1));
    tempToOR.resize(numVectors * (delta-1));
    for (int j=0; j<numVectors; j++) {
        for (int i = 0; i < delta - 1; i++) {
            sharesVector[j*(delta - 1) + i].resize(i + 2);
            memcpy((byte *) sharesVector[j*(delta - 1) + i].data(), (byte *) &a[j*delta], (i + 2) * field->getElementSizeInBytes());
            totalSize += i + 2;
        }
        b[j*delta] = a[j*delta];
    }
    unboundedOR(sharesVector, totalSize, tempToOR, 0);
    for (int j=0; j<numVectors; j++) {
        for (int i = 0; i < delta - 1; i++) {
            b[j*delta + i + 1] = tempToOR[j*(delta-1) + i];
        }
    }
//    b[0] = a[0];
//    for (int i=1; i<delta; i++){
//        b[i] = unboundedOR(a.data(), i+1);
//    }

//    temp.resize(b.size());
//    openShare(b,temp,T);
//    cout<<"b:"<<endl;
//    for (int i=0; i<temp.size(); i++){
//        cout<<temp[i]<<" ";
//    }
//    cout<<endl;

    //calculate si = all zeros including the group with the first "1" bit and then all ones.
    for (int j=0; j<numVectors; j++) {
        for (int i = 0; i < lambda; i++) {
            s[j*lambda + i] = y[j*lambda + i] - f[j*lambda + i];

        }
    }

//    temp.resize(s.size());
//    openShare(s,temp,T);
//    cout<<"s:"<<endl;
//    for (int i=0; i<temp.size(); i++){
//        cout<<temp[i]<<" ";
//    }
//    cout<<endl;

    vector<FieldType> fToMultiply(numVectors*delta*lambda);
    vector<FieldType> bToMultiply(numVectors*delta*lambda);
    //calculate the real b values:
    for (int k = 0; k<numVectors; k++) {
        for (int i = 0; i < lambda; i++) {
            for (int j = 0; j < delta; j++) {
                fToMultiply[k*delta*lambda + i * delta + j] = f[k*lambda + i];
                bToMultiply[k*delta*lambda + i * delta + j] = b[k*delta + j];
            }
        }
    }

    DNHonestMultiplication(fToMultiply.data(), bToMultiply.data(), fToMultiply, numVectors*delta*lambda);

//        cout<<"fToMultiply:"<<endl;
//        for (int i=0; i<delta; i++){
//            cout<<fToMultiply[i]<<" ";
//        }
//        cout<<endl;

    for (int k = 0; k<numVectors; k++) {
        for (int i = 0; i < lambda; i++) {
            if (i < lambda - 1) {
                for (int j = 0; j < delta; j++) {

                    bToFill[k][i * delta + j] = fToMultiply[k*lambda*delta  + i * delta + j] + s[k*lambda + i];
                }
            } else { //last group
                int remain = shares[k].size() % delta;
                for (int j = 0; j < remain; j++) {
                    bToFill[k][i * delta + j] = fToMultiply[k*lambda*delta  + i * delta + j] + s[k*lambda + i];
                }
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
//
//template <class FieldType>
//void ProtocolParty<FieldType>::prefixOR(FieldType* shares, vector<FieldType> & bToFill) {
//    int size = bToFill.size();
////    cout<<"original size = "<<size<<endl;
//
//    //calculate xi = unbounded or of each delta bits
//    vector<FieldType> x(lambda);
//    vector<FieldType> y(lambda);
//    vector<FieldType> f(lambda);
//    vector<FieldType> s(lambda);
//
//    vector<vector<FieldType>> sharesVector(lambda);
//
//    int totalSize;
//    //copy the first lambda-1 groups
//    for (int i=0; i<lambda-1; i++){
//        sharesVector[i].resize(delta);
//        memcpy((byte*)sharesVector[i].data(), (byte*)&shares[i*delta], delta*field->getElementSizeInBytes());
//    }
//    //copy the last gropu. Notice that it can contain less than delta elements
//    int remain = size % delta;
//    if (remain == 0){
//        sharesVector[lambda-1].resize(delta);
//        memcpy((byte*)sharesVector[lambda-1].data(), (byte*)&shares[(lambda-1)*delta], delta*field->getElementSizeInBytes());
//        totalSize = lambda*delta;
//
//    } else {
//        sharesVector[lambda-1].resize(remain);
//        memcpy((byte*)sharesVector[lambda-1].data(), (byte*)&shares[(lambda-1)*delta], remain*field->getElementSizeInBytes());
//        totalSize = (lambda-1)*delta + remain;
//    }
//
//    vector<FieldType> temp;
////    cout<<"sharesAlligned:"<<endl;
////    for (int i=0; i<lambda; i++){
////        temp.resize(sharesVector[i].size());
////        openShare(sharesVector[i],temp,T);
////        for (int i=0; i<temp.size(); i++){
////            cout<<temp[i]<<" ";
////        }
////    }
////    cout<<endl;
////cout<<"before unbounded"<<endl;
////    cout<<"totalSize = "<<totalSize<<endl;
//    unboundedOR(sharesVector, totalSize, x, 0);
////    for (int i=0; i<lambda; i++){
////        x[i] = unboundedOR(sharesAlligned + i*delta, delta);
////    }
//
////    temp.resize(x.size());
////    openShare(x,temp,T);
////    cout<<"x:"<<endl;
////    for (int i=0; i<temp.size(); i++){
////        cout<<temp[i]<<" ";
////    }
////    cout<<endl;
//
//    //calculate yi = prefix or of x
//    totalSize = 0;
//    sharesVector.resize(lambda-1);
//    for (int i=0; i<lambda-1; i++){
//        sharesVector[i].resize(i+2);
//        memcpy((byte*)sharesVector[i].data(), (byte*)x.data(), (i+2)*field->getElementSizeInBytes());
//        totalSize += i+2;
//    }
//    y[0] = x[0];
//    unboundedOR(sharesVector, totalSize, y, 1);
////    for (int i=1; i<lambda; i++){
////        y[i] = unboundedOR(x.data(), i+1);
////    }
//
////    openShare(y,temp,T);
////    cout<<"y:"<<endl;
////    for (int i=0; i<temp.size(); i++){
////        cout<<temp[i]<<" ";
////    }
////    cout<<endl;
//
//    //calculate fi = find the group with the first "1" bit
//    f[0] = y[0];
//    for (int i=1; i<lambda; i++){
//        f[i] = y[i] - y[i-1];
//    }
//
//    //find the real b values for the group which contains the first "1" bit.
//    //Note that the group index is unknown
//    vector<FieldType> a(delta, *field->GetZero());
//    vector<FieldType> aToMultiply(lambda);
//
//    //TODO can be optimized!!
//    for (int j=0; j<delta; j++){
////        for (int i=0; i<lambda; i++) {
////            aToMultiply[i] = sharesAlligned[i * delta + j];
////        }
//
//        for (int i=0; i<lambda-1; i++) {
//            a[j] += f[i]*shares[i * delta + j];
//
//        }
//
//        if (j<remain){
//            a[j] += f[lambda-1]*shares[(lambda-1) * delta + j];
//        }
//    }
//
//    DNHonestMultiplication(a.data(), a.data(), a, delta, true);
//
//    //calculate the real b values for the group which contains the first "1" bit.
//    //Note that the group index is unknown
//    vector<FieldType> b(delta);
//
//    totalSize = 0;
//    sharesVector.resize(delta-1);
//    for (int i=0; i<delta-1; i++){
//        sharesVector[i].resize(i+2);
//        memcpy((byte*)sharesVector[i].data(), (byte*)a.data(), (i+2)*field->getElementSizeInBytes());
//        totalSize += i+2;
//    }
//    b[0] = a[0];
//    unboundedOR(sharesVector, totalSize, b, 1);
////    b[0] = a[0];
////    for (int i=1; i<delta; i++){
////        b[i] = unboundedOR(a.data(), i+1);
////    }
//
//    //calculate si = all zeros including the group with the first "1" bit and then all ones.
//    for (int i=0; i<lambda; i++){
//        s[i] = y[i] - f[i];
//
//    }
//
//    vector<FieldType> fToMultiply(delta*lambda);
//    vector<FieldType> bToMultiply(delta*lambda);
//    //calculate the real b values:
//    for (int i=0; i<lambda; i++) {
//        for (int j = 0; j < delta; j++) {
//            fToMultiply[i*delta + j] = f[i];
//            bToMultiply[i*delta + j] = b[j];
//        }
//    }
//
//    DNHonestMultiplication(fToMultiply.data(), bToMultiply.data(), fToMultiply, delta*lambda);
//
////        cout<<"fToMultiply:"<<endl;
////        for (int i=0; i<delta; i++){
////            cout<<fToMultiply[i]<<" ";
////        }
////        cout<<endl;
//
//    for (int i=0; i<lambda; i++) {
//        if (i<lambda-1) {
//            for (int j = 0; j < delta; j++) {
//
//                bToFill[i * delta + j] = fToMultiply[i * delta + j] + s[i];
//            }
//        } else { //last group
//            for (int j = 0; j < remain; j++) {
//                bToFill[i * delta + j] = fToMultiply[i * delta + j] + s[i];
//            }
//        }
//    }
//
////    openShare(x,x,T);
////    cout<<"x:"<<endl;
////    for (int i=0; i<lambda; i++){
////        cout<<x[i]<<" ";
////    }
////    cout<<endl;
////    openShare(y,y,T);
////    cout<<"y:"<<endl;
////    for (int i=0; i<lambda; i++){
////        cout<<y[i]<<" ";
////    }
////    cout<<endl;
////    openShare(f,f,T);
////    cout<<"f:"<<endl;
////    for (int i=0; i<lambda; i++){
////        cout<<f[i]<<" ";
////    }
////    cout<<endl;
////
////    openShare(b,b,T);
////    cout<<"b_i0:"<<endl;
////    for (int i=0; i<delta; i++){
////        cout<<b[i]<<" ";
////    }
////    cout<<endl;
////
////    openShare(s,s,T);
////    cout<<"s:"<<endl;
////    for (int i=0; i<lambda; i++){
////        cout<<s[i]<<" ";
////    }
////    cout<<endl;
//
//
//}


template <class FieldType>
void ProtocolParty<FieldType>::bitwiseLessThan(vector<vector<FieldType>> & a, vector<vector<FieldType>> & b, vector<FieldType> & result) {
    int numVectors = a.size();
    int size = a[0].size();
//    vector<FieldType> mults(size);
//    DNHonestMultiplication(a, b, mults, size);
    FieldType two(2);

    //calculate [ci] = [ai]^[bi]
    vector<vector<FieldType>> c(numVectors, vector<FieldType>(size));
    for (int j=0; j<numVectors; j++) {
        for (int i = 0; i < size; i++) {
            c[j][i] = a[j][i] + b[j][i];
            c[j][i] = c[j][i] - (two * a[j][i] * b[j][i]);
        }
    }

//    cout<<" in bitwise"<<endl;
//    vector<FieldType> temp(size);
//    for (int j=0; j<numVectors;j++) {
//        openShare(c[j], temp, T);
//        cout << "c:" << endl;
//        for (int i = 0; i < size; i++) {
//            cout << temp[i] << " ";
//        }
//        cout << endl;
//    }

    //compute d[i] = prefix OR of ci
    vector<vector<FieldType>> d(numVectors, vector<FieldType>(size));
    auto start = high_resolution_clock::now();
    prefixOR(c, d);
    auto end = high_resolution_clock::now();
    debugDuration += duration_cast<microseconds>(end- start).count();

//    for (int j=0; j<numVectors;j++) {
//        openShare(d[j], temp, T);
//        cout << "d:" << endl;
//        for (int i = 0; i < size; i++) {
//            cout << temp[i] << " ";
//        }
//        cout << endl;
//    }

    //compute e = [di]-[di+1]
    vector<FieldType> e(numVectors*size);
    for (int j=0; j<numVectors; j++) {
        e[j*size] = d[j][0];
        for (int i = 1; i < size; i++) {
            e[j*size + i] = d[j][i] - d[j][i - 1];
        }
    }
//    openShare(e,temp,T);
//    cout<<"e:"<<endl;
//    for (int i=0; i<size; i++){
//        cout<<temp[i]<<" ";
//    }
//    cout<<endl;

    vector<FieldType> res(numVectors, *field->GetZero());
    for (int j=0; j<numVectors; j++) {
        for (int i = 0; i < size; i++) {
            res[j] += e[j*size + i] * b[j][i];
//        res += mults[i];
        }
    }

//    vector<FieldType> mults(size);
    DNHonestMultiplication(res.data(), res.data(), result, numVectors, true);
}

template <class FieldType>
void ProtocolParty<FieldType>::geneateRandomWithBits(int size, vector<FieldType> & randomToFill, vector<vector<FieldType>> & bitsToFill, int numRandoms) {

    int fullSize = numRandoms*1.1;
    randomToFill.resize(fullSize);
    bitsToFill.resize(fullSize);

    vector<FieldType> toCheck(fullSize);
    vector<FieldType> rToCheck(fullSize);
    vector<FieldType> checkResult(fullSize);
    FieldType p(FieldType::p);

    for (int j=0; j<fullSize; j++) {
        bitsToFill[j].resize(size);
        randomBit(bitsToFill[j]);

        randomToFill[j] = *field->GetZero();
        for (int i = 0; i < size; i++) {
            randomToFill[j] += bitsToFill[j][i] * twoSquares[twoSquares.size() - size + i];
        }

        toCheck[j] = p - randomToFill[j];
        toCheck[j] = toCheck[j] * randomSharesArray[randomSharesOffset++];
    }

    openShare(toCheck, toCheck, 2*T);
    vector<int> wrongValuesIndices;
    vector<int> rightValuesIndicesInSpare;
    for (int j=0; j<fullSize; j++) {
        if (j<numRandoms) {
            if (toCheck[j] == *field->GetZero()) {
                wrongValuesIndices.push_back(j);
            }
        } else {
            if (toCheck[j] != *field->GetZero()) {
                rightValuesIndicesInSpare.push_back(j);
            }
        }
    }

    int counter = 0;
    for (int j=0; j<wrongValuesIndices.size(); j++) {
        randomToFill[wrongValuesIndices[j]] = randomToFill[rightValuesIndicesInSpare[counter]];
        memcpy((byte*) bitsToFill[wrongValuesIndices[j]].data(), (byte*) bitsToFill[rightValuesIndicesInSpare[counter++]].data(), size*field->getElementSizeInBytes());
    }


}

template <class FieldType>
void ProtocolParty<FieldType>::LSB(vector<FieldType> & x, int bitsSize, vector<FieldType> & lsb) {

    int numElements = x.size();
//    vector<vector<FieldType>> rBits = bitsToLSB[randomWithBitsOffset];
//    vector<FieldType> r = randomsToLSB[randomWithBitsOffset++];

//    geneateRandomWithBits(bitsSize, r, rBits);
    vector<FieldType> temp(bitsSize);
//for (int j=0; j<numElements;j++) {
//
////    openShare(bitsToLSB[j], temp, T);
////    cout << "bits of random:" << endl;
////    for (int i = 0; i < bitsSize; i++) {
////        cout << temp[i];
////    }
////    cout << endl;
//    temp[0] = randomsToLSB[j];
//    openShare(temp, temp, T);
//    cout << "r:" << temp[0] << endl;
//}


    //TODO currently we use just the 3 generates randoms...
    vector<FieldType> c(numElements);
    vector<FieldType> cOpened(numElements);
    for (int i=0; i<numElements; i++) {
        c[i] = x[i] + randomsToLSB[i];
    }
    openShare(c, cOpened, T);
//    cout<<"c:"<<cOpened[0]<<endl;
//    cout<<"c:"<<cOpened[1]<<endl;
//    cout<<"c:"<<cOpened[2]<<endl;

    vector<vector<FieldType>> cBits(numElements, vector<FieldType>(bitsSize));
    getBits(cOpened, cBits, bitsSize);

//    for (int j=0; j<numElements; j++) {
//    openShare(cBits[j], temp, T);
//    cout<<"c bits:"<<endl;
//
//        for (int i = 0; i < bitsSize; i++) {
//            cout << temp[i];
//        }
//        cout << endl;
//
//    }

    FieldType two(2);
    vector<FieldType> compare(numElements);
    bitwiseLessThan(cBits, bitsToLSB, compare);

//    openShare(compare, temp, T);

//    for (int j=0; j<numElements; j++) {
//        cout<<"c<r ? "<< temp[j]<<endl;
//
//    }
//    temp[0] = compare;
//    openShare(temp, temp, T);
//     cout<<"c<r? "<<temp[0]<<endl;
    vector<FieldType> left(2*numElements);
    vector<FieldType> right(2*numElements);
    vector<FieldType> mults(2*numElements);

    for (int i=0; i<numElements; i++) {
        FieldType xorBits0 = cBits[i][bitsSize - 1] + bitsToLSB[i][bitsSize - 1] - two * cBits[i][bitsSize - 1] * bitsToLSB[i][bitsSize - 1];
//    temp[0] = xorBits0;
//    openShare(temp, temp, T);
//    cout<<"c0 ^ r0 = "<<temp[0]<<endl;
        left[i*2] = compare[i];
        right[i*2] = *field->GetOne() - xorBits0;
        left[i*2 + 1] = *field->GetOne() - compare[i];
        right[i*2 + 1] = xorBits0;
    }

    DNHonestMultiplication(left.data(), right.data(), mults, 2*numElements);

    for (int i=0; i<numElements; i++) {
        lsb[i] = mults[i*2] + mults[i*2 + 1];
    }
//    temp[0] = res;
//    openShare(temp, temp, T);
//    cout<<"res = "<<temp[0]<<endl;
//    return res;

}

template <class FieldType>
FieldType ProtocolParty<FieldType>::compare(FieldType & a, FieldType & b, int bitSize){


    FieldType two(2);
    vector<FieldType> toLSB(3);
    vector<FieldType> lsb(3);
    toLSB[0] = a*two;
    toLSB[1] = b*two;
    toLSB[2] = (a-b)*two;
    LSB(toLSB, bitSize, lsb);
    //w = [a<p/2]
    vector<FieldType> w(1);
    w[0] = *field->GetOne() - lsb[0];

    //x = [b<p/2]
    vector<FieldType> x(1);
    x[0] = *field->GetOne() - lsb[1];

    //y = [(a-b)<p/2]
    vector<FieldType> y(1);
    y[0] = *field->GetOne() - lsb[2];

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
void ProtocolParty<FieldType>::getBits(vector<FieldType> & x, vector<vector<FieldType>> & bits, int size) {
    int numElements = x.size();
    for (int j=0; j<numElements; j++) {
        for (int i = size - 1; i >= 0; i--) {
            bits[j][i] = x[j].elem & 1;
            x[j].elem = x[j].elem >> 1;
        }
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
    int max = delta > lambda ? delta : lambda;
    vector<FieldType> x(max+1);
    vector<FieldType> y(max+1);

    polynomials.resize(max-1);

    for (int i = 0; i < max+1; i++) {
        x[i] = i + 1;
        y[i] = (i == 0) ? 0 : 1;
    }
    for (int j=2; j<=max; j++) {
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

    numCompares = 100;
    debugDuration = 0;

}


template <class FieldType>
bool ProtocolParty<FieldType>::preparationPhase() {

    int numOfRandomShares = 2*times*numCompares*3*(2*(lambda*delta + (delta*(delta-1))/2 - 1 + (lambda*(lambda-1))/2 - 1) + field->getElementSizeInBits());
    randomSharesOffset = 0;
    randomSharesArray.resize(numOfRandomShares);

    //generate enough random shares for the AES key
//    generateRandomShares(numOfRandomShares, randomSharesArray);
    RandomSharingForInputs(numOfRandomShares, randomSharesArray);

    //run offline for all the future multiplications including the multiplication of the protocol

    int numOfRandomSharesForMults = 2*times*numCompares*3*(lambda*delta + (delta*(delta-1))/2 - 1 + (lambda*(lambda-1))/2 - 1 + 2 + 1 + 1 + lambda*delta);
    randomTAnd2TSharesOffset = 0;
    generateRandom2TAndTShares(numOfRandomShares,randomTAnd2TShares);

    randomWithBitsOffset = 0;
    geneateRandomWithBits(field->getElementSizeInBits(), randomsToLSB, bitsToLSB, 3*numCompares);
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
void ProtocolParty<FieldType>::DNHonestMultiplication(FieldType *a, FieldType *b, vector<FieldType> &cToFill, int numOfTupples, bool sop) {

    int fieldByteSize = field->getElementSizeInBytes();
    vector<FieldType> xyMinusRShares(numOfTupples);//hold both in the same vector to send in one batch
    vector<byte> xyMinusRSharesBytes(numOfTupples *fieldByteSize);//hold both in the same vector to send in one batch

    vector<FieldType> xyMinusR;//hold both in the same vector to send in one batch
    vector<byte> xyMinusRBytes;

    vector<vector<FieldType>> recBufsElements(N);
    vector<vector<FieldType>> sendBufsElements(N);

    if (!sop) {
        //generate the shares for x+a and y+b. do it in the same array to send once
        for (int k = 0; k < numOfTupples; k++)//go over only the logit gates
        {
            //compute the share of xy-r
            xyMinusRShares[k] = a[k] * b[k] - randomTAnd2TShares[randomTAnd2TSharesOffset + 2 * k + 1];

        }
    } else { //a is sum of products. No need to multiply
        //generate the shares for x+a and y+b. do it in the same array to send once
        for (int k = 0; k < numOfTupples; k++)//go over only the logit gates
        {
            //compute the share of xy-r
            xyMinusRShares[k] = a[k] - randomTAnd2TShares[randomTAnd2TSharesOffset + 2 * k + 1];

        }
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
void ProtocolParty<FieldType>::openShare(vector<FieldType> &shares, vector<FieldType> &secrets,
                                                  int d){

    int numOfRandomShares = shares.size();
    //vector<vector<byte>> sendBufsBytes(N);

    vector<vector<FieldType>> sendBufsElements(N);
    vector<vector<FieldType>> recBufsElements(N,vector<FieldType>(numOfRandomShares));

    vector<FieldType> x1(N);


    //copy the same data for all parties
    for(int i=0; i<N; i++){

        sendBufsElements[i] = shares;
    }

    //call the round function to send the shares to all the users and get the other parties share
    roundFunctionSyncElements(sendBufsElements, recBufsElements,12);

    //reconstruct each set of shares to get the secret

    for(int k=0; k<numOfRandomShares; k++){

        //get the set of shares for each element
        for(int i=0; i < N; i++) {

            x1[i] = recBufsElements[i][k];
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
FieldType ProtocolParty<FieldType>::evalPolynomial(const vector<FieldType>& coeff, FieldType* a)
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
