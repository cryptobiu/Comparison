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

#define flag_print true
#define flag_print_timings true

template <class FieldType>
class ProtocolParty : public MPCProtocol, public HonestMajority {

private:
    int N, T, partyId;//number of servers

    int times; //number of times to run the run function
    int iteration; //number of the current iteration

    boost::asio::io_service io_service;
    vector<PrgFromOpenSSLAES> prgs;

    TemplateField<FieldType> *field;
    vector<shared_ptr<ProtocolPartyData>> parties;
    Measurement* timer;

    vector<FieldType> alpha; // N distinct non-zero field elements
    vector<FieldType> beta;
    VDM<FieldType> matrix_vand;
    VDMTranspose<FieldType> matrix_vand_transpose;
    HIM<FieldType> matrix_for_interpolate;
    vector<FieldType> y_for_interpolate;

    vector<FieldType> randomTAnd2TShares;
    int randomTAnd2TSharesOffset = 0;
    vector<FieldType> randomSharesArray;
    int randomSharesOffset = 0;


    void initializationPhase();
    bool preparationPhase();

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
//
//
//    shiftbyOne.resize(securityParamter);
//    for(int i=0; i<securityParamter; i++){
//        shiftbyOne[i] = 1 << i;
//    }
//
//

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

    cout<<"*********** check inverse *************"<<endl;
    vector<FieldType> element(10);
    vector<FieldType> inverseArr(10);
    element[0] = randomSharesArray[randomSharesOffset++];
    cout<<"share of element is " << element[0]<<endl;

    vector<FieldType> openedElement(10);
    openShare(element, openedElement, T);
    cout<<"opened element is "<<openedElement[0]<<endl;

    inverse(element.data(), inverseArr);
    cout<<"element inverse is " << inverseArr[0]<<endl;


    vector<FieldType> openedInverse(10);
    openShare(inverseArr, openedInverse, T);
    cout<<"opened element inverse is "<<openedInverse[0]<<endl;

    for (int i=0; i<10; i++) {
        auto res = openedElement[i] * openedInverse[i];

        cout << "the output is " << res << endl;
    }

    cout<<"*********** check RANDOM BIT *************"<<endl;
    vector<FieldType> bits(10);
    vector<FieldType> openBits(10);
    randomBit(bits);

    cout<<"shares of random bits:"<<endl;
    for (int i=0; i<10; i++){
        cout<<bits[i]<<" ";
    }
    cout<<endl;
    openShare(bits, openBits, T);
    cout<<"actual random bits:"<<endl;
    for (int i=0; i<10; i++){
        cout<<openBits[i]<<" ";
    }
    cout<<endl;

//    auto t1 = high_resolution_clock::now();
//    timer->startSubTask("inputPhase", iteration);
//    inputPhase();
//    timer->endSubTask("inputPhase", iteration);
//    auto t2 = high_resolution_clock::now();
//    auto duration = duration_cast<milliseconds>(t2-t1).count();
//
//    if(flag_print_timings) {
//        cout << "time in milliseconds inputPhase: " << duration << endl;
//    }



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
//    matrix_him.allocate(N,N,field);
    matrix_vand.allocate(N,N,field);
    matrix_vand_transpose.allocate(N,N,field);
//    m.allocate(T, N-T,field);

    // Compute Vandermonde matrix VDM[i,k] = alpha[i]^k
    matrix_vand.InitVDM();
    matrix_vand_transpose.InitVDMTranspose();

    // Prepare an N-by-N hyper-invertible matrix
//    matrix_him.InitHIM();

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


}


template <class FieldType>
bool ProtocolParty<FieldType>::preparationPhase() {

    int keysize = 16/field->getElementSizeInBytes() + 1;

    int numOfRandomShares = 10*keysize + 1;
    randomSharesArray.resize(numOfRandomShares);

    //generate enough random shares for the AES key
    generateRandomShares(numOfRandomShares, randomSharesArray);


    //run offline for all the future multiplications including the multiplication of the protocol

    randomTAnd2TSharesOffset = 0;
    generateRandom2TAndTShares(numParties*100,randomTAnd2TShares);


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
void ProtocolParty<FieldType>::inverse(FieldType *a, vector<FieldType> &bToFill) {

    int numOfElements = bToFill.size();
    vector<FieldType> raShares(numOfElements);//hold both in the same vector to send in one batch
    vector<FieldType> openedRAShares(numOfElements);
    FieldType one(1);


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
        bToFill[k] = (one / openedRAShares[k]) * randomSharesArray[randomSharesOffset + k];
//        cout<<"inverse of a share = "<< bToFill[k]<<endl;
    }

    randomSharesOffset += numOfElements;

}


template <class FieldType>
void ProtocolParty<FieldType>::randomBit(vector<FieldType> &bToFill) {

    int numOfElements = bToFill.size();
    cout<<"numOfElements = "<<numOfElements<<endl;
    vector<FieldType> rSquareShares(numOfElements);//hold both in the same vector to send in one batch
    vector<FieldType> openedRSquareShares(numOfElements);
    FieldType one(1);
    FieldType twoInv(2);
    twoInv = one / twoInv;

    FieldType r;
    cout<<"randomSharesOffset = "<<randomSharesOffset<<endl;
    //generate the shares for x+a and y+b. do it in the same array to send once
    for (int k = 0; k < numOfElements; k++)//go over only the logit gates
    {
        //compute the share of r*a
        r = randomSharesArray[randomSharesOffset + k];
        rSquareShares[k] = r*r;
    }

    openShare(rSquareShares, openedRSquareShares, 2*T);

    vector<FieldType> roots(numOfElements);
    for (int k = 0; k < numOfElements; k++)//go over only the logit gates
    {
        roots[k] = openedRSquareShares[k].sqrt();

    }

    vector<FieldType> rootsInverse(numOfElements);
    inverse(roots.data(), rootsInverse);

    for (int k = 0; k < numOfElements; k++)
    {
        bToFill[k] = roots[k] * rootsInverse[k];
        bToFill[k] = (bToFill[k] + one) * twoInv;
    }

    randomSharesOffset += numOfElements;

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


#endif //COMPARISON_PROTOCOLPARTY_H
