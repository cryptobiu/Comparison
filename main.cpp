#include <iostream>
#include "ProtocolParty.h"
#include "SmallFields.h"

int main(int argc, char* argv[])
{
    CmdParser parser;
    auto parameters = parser.parseArguments("", argc, argv);
    int times = stoi(parser.getValueByKey(parameters, "internalIterationsNumber"));


    string fieldType = parser.getValueByKey(parameters, "fieldType");
    cout<<"fieldType = "<<fieldType<<endl;

    if(fieldType.compare("ZpMersenne31") == 0)
    {
        ProtocolParty<ZpMersenneIntElement> protocol(argc, argv);
        auto t1 = high_resolution_clock::now();
        protocol.run();

        auto t2 = high_resolution_clock::now();

        auto duration = duration_cast<milliseconds>(t2-t1).count();
        cout << "time in milliseconds for " << times << " runs: " << duration << endl;
        cout << "end main" << '\n';

    }
    else if(fieldType.compare("ZpMersenne13") == 0)
    {

        ProtocolParty<ZpMersenne13> protocol(argc, argv);
        auto t1 = high_resolution_clock::now();
        protocol.run();
        auto t2 = high_resolution_clock::now();

        auto duration = duration_cast<milliseconds>(t2-t1).count();
        cout << "time in milliseconds for " << times << " runs: " << duration << endl;
        cout << "end main" << '\n';

    }

    else if(fieldType.compare("Zp16BitPrime") == 0)
    {

        ProtocolParty<Zp16BitPrime> protocol(argc, argv);
        auto t1 = high_resolution_clock::now();
        protocol.run();
        auto t2 = high_resolution_clock::now();

        auto duration = duration_cast<milliseconds>(t2-t1).count();
        cout << "time in milliseconds for " << times << " runs: " << duration << endl;
        cout << "end main" << '\n';

    }
    return 0;
}
