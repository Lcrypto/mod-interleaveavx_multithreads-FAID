#include "CChannel.h"
#include "CLDPC.h"
#include "CModulate.h"
#include "CSimulate.h"
#include "CTool.h"
#include <fstream>
#include <iomanip>
#include <iostream>
#include <stdio.h>
#include <sys/time.h> // Used for calculation Time
#include <time.h>
#include <thread>
using namespace std;
int collectflag;
int MAX_THREADS;

int main(void)
{
    int MAX_THREADS = std::thread::hardware_concurrency();
    cout << "sizeof(CLDPC)=" << sizeof(CLDPC) << endl;
    cout << "sizeof(CMoudlate)=" << sizeof(CModulate) << endl;
    cout << "sizeof(CChannel)=" << sizeof(CChannel) << endl;
    cout << "sizeof(CSimulate)=" << sizeof(CSimulate) << endl;
    cout << "sizeof(Statistic)=" << sizeof(Statistic) << endl;
    int i;
    timeval t_start, t_end;
    Parameter_Simulation p_simulation;
    ReadProfile(&p_simulation);

    /*MAX_THREADS threads initialize*/
    CSimulate simulate[MAX_THREADS];
    for (i = 0; i < MAX_THREADS; ++i) {
        simulate[i].Initial(p_simulation, i);
    }

    float snr = p_simulation.snr_start;
    float snr_end = p_simulation.snr_end;

    unsigned long TestFrame = 0;
    unsigned long ErrorFrame = 0;
    unsigned long ErrorBits = 0;
    unsigned long LT3ErrBitFrame = 0;
    unsigned long ModErrorBits = 0;
    unsigned long ModErrorSymbol = 0;
    unsigned long ModErrorFrame = 0;
    double BER = 0;
    double ModBER = 0;
    double FER = 0;
    double ModSER = 0;
    double ModFER = 0;

    ofstream fout; //("Result.txt", std::ios::app);
    ofstream tout;
    ofstream eout;
    ofstream nout;
    ofstream dout;
    ofstream demodout;
    ofstream iterOut;

    collectflag = 0;
    fout.open("Result.txt", std::ios::app);
    if (!fout.is_open()) {
        cerr << "Cannot open Result.txt\n";
        getchar();
        exit(EXIT_FAILURE);
    }

    tout.open("Temp.txt", std::ios::app);
    if (!tout.is_open()) {
        cerr << "Cannot open Temp.txt\n";
        getchar();
        exit(EXIT_FAILURE);
    }

    demodout.open("demod.txt", std::ios::app);
    if (!demodout.is_open()) {
        cerr << "Cannot open Temp.txt\n";
        getchar();
        exit(EXIT_FAILURE);
    }
    demodout << setw(5) << " Eb/N0" << '\t' << setw(20) << " ModFER" << '\t' << setw(20) << "ModBER" << '\t' << setw(20)
             << "ModSER" << '\t' << endl;
    demodout.close();
    tout.open("Temp.txt", std::ios::app);
    /* Basic information of Temp.txt file will be lost after refresh*/
    tout << "DATE:" << __DATE__ << endl;
    tout << "Time" << __TIME__ << endl;
    tout << "MaxItertion:" << p_simulation.Max_Iteration << endl;
    tout << "Modulation Type" << p_simulation.mod_type << endl;
    tout << "InterLeave ModType" << p_simulation.interleavemod_type << endl;
    tout << "scale=" << p_simulation.scale << endl;
    tout << "codeFile:" << p_simulation.fileName << endl;
    tout << "Punctue Number: " << _PunctureBits << " Shorten Bits" << _ShortenBits
         << " RATE: " << simulate[0].ldpc->m_Rate << endl;
    tout << setw(5) << "Eb_N0" << setw(20) << "TestFrame" << setw(15) << "ErrorFrame" << setw(20) << "ErrorBits"
         << setw(20) << "FER" << setw(20) << "BER" << setw(15) << "LT3ErrBitFrame" << setw(15) << endl;
    tout.close();

    /* Basic information of Result.txt file will be lost after later refresh*/
    fout << endl;
    fout << "**********************************************************************************************************"
            "**********************************"
         << endl;
    fout << "DATE:" << __DATE__ << endl;
    fout << "Time" << __TIME__ << endl;
    fout << "codeFile:" << p_simulation.fileName << endl;
    fout << "DecodeMethod:" << p_simulation.decode_method << endl;
    fout << "MaxItertion:" << p_simulation.Max_Iteration << endl;
    fout << "Modulation Type" << p_simulation.mod_type << endl;
    fout << "InterLeave ModType" << p_simulation.interleavemod_type << endl;
    fout << "scale=" << p_simulation.scale << endl;
    fout << "factor_1=" << p_simulation.Factor_1 << endl;
    fout << "factor_2=" << p_simulation.Factor_2 << endl;

    fout << "Punctue Number: " << _PunctureBits << " Shorten Bits" << _ShortenBits
         << " RATE: " << simulate[0].ldpc->m_Rate << endl;
    fout << setw(5) << "Eb_N0" << '\t' << setw(20) << "TestFrame" << '\t' << setw(15) << "ErrorFrame" << '\t'
         << setw(20) << "ErrorBits" << '\t' << setw(20) << "FER" << '\t' << setw(20) << "BER" << '\t' << setw(15)
         << "LT3ErrBitFrame" << '\t' << setw(15) << "Time(s)" << '\t' << endl;
    cout << "DATE:" << __DATE__ << endl;
    cout << "Time" << __TIME__ << endl;
    cout << "codeFile:" << p_simulation.fileName << endl;
    cout << "DecodeMethod:" << p_simulation.decode_method << endl;
    cout << "MaxItertion:" << p_simulation.Max_Iteration << endl;
    cout << "Modulation Type" << p_simulation.mod_type << endl;
    cout << "InterLeave ModType" << p_simulation.interleavemod_type << endl;
    cout << "scale=" << p_simulation.scale << endl;
    cout << "factor_1=" << p_simulation.Factor_1 << endl;
    cout << "factor_2=" << p_simulation.Factor_2 << endl;
    cout << "Punctue Number: " << _PunctureBits << " Shorten Bits" << _ShortenBits
         << " RATE: " << simulate[0].ldpc->m_Rate << endl;
    cout << setw(5) << "Eb_N0" << setw(20) << "TestFrame" << setw(15) << "ErrorFrame" << setw(20) << "ErrorBits"
         << setw(20) << "FER" << setw(20) << "BER" << setw(15) << "LT3ErrBitFrame" << setw(15) << "Time(s)" << endl;
    fout.close();

    for (snr = p_simulation.snr_start; snr < snr_end; snr += p_simulation.snr_pass) {
        TestFrame = 0;
        ErrorFrame = 0;
        ErrorBits = 0;
        ModErrorBits = 0;
        ModErrorFrame = 0;
        ModErrorSymbol = 0;
        LT3ErrBitFrame = 0;
        demodout.open("demod.txt", std::ios::app);
        fout.open("Result.txt", std::ios::app); //Use ios::app|ios::out to create a file if it does not exist, or append it to the end of the file if it does exist.
        eout.open("errorindex.txt", std::ios::app);
        nout.open("errorfloat.txt", std::ios::app);
        dout.open("errordecode.txt", std::ios::app);
        iterOut.open("iterCount.txt", std::ios::app);
        iterOut << "Eb/N0: " << setw(5) << snr << "scale=" << p_simulation.scale << endl;
        eout << "Eb/N0: " << setw(5) << snr << "scale=" << p_simulation.scale << endl;
        nout << "Eb/N0: " << setw(5) << snr << "scale=" << p_simulation.scale << endl;
        dout << "Eb/N0: " << setw(5) << snr << "scale=" << p_simulation.scale << endl;
        eout.close();
        nout.close();
        dout.close();
        iterOut.close();
        //Initialization noise
        for (i = 0; i < MAX_THREADS; ++i) {
            simulate[i].Configure(snr, p_simulation.decode_method);
        }
        gettimeofday(&t_start, NULL);

        while (TestFrame < 1000 || ErrorFrame < 20) {

            for (i = 0; i < MAX_THREADS; ++i) {
                simulate[i].Start();
            }

            for (i = 0; i < MAX_THREADS; ++i) {
                simulate[i].End();
            }

            for (i = 0; i < MAX_THREADS; ++i) {
                TestFrame += simulate[i].TestFrame;
                ErrorFrame += simulate[i].ErrorFrame;
                ErrorBits += simulate[i].ErrorBits;
                LT3ErrBitFrame += simulate[i].LT3ErrBitFrame;
                ModErrorBits += simulate[i].ModErrorBits;
                ModErrorFrame += simulate[i].ModErrorFrame;
                ModErrorSymbol += simulate[i].ModErrorSymbol;
            }
            ModBER = (double)(ModErrorBits) / (TestFrame * (NmoinsK - _ShortenBits));
            ModSER = (double)(ModErrorSymbol) / (TestFrame * (NmoinsK - _ShortenBits) / p_simulation.mod_type);
            ModFER = (double)(ModErrorFrame) / TestFrame;
            BER = (double)(ErrorBits > 0 ? ErrorBits : 1)
                / (TestFrame * (NmoinsK - _ShortenBits)); // Assume one is wrong, simulate BER
            FER = (double)(ErrorFrame > 0 ? ErrorFrame : 1) / TestFrame;

            if (FER < 1E-5) {
                collectflag = 1;
            }

            tout.open("Temp.txt", std::ios::out); // ios::out file is opened in output mode (memory data is output to the file)

            tout << setw(5) << snr << '\t' << setw(20) << TestFrame << '\t' << setw(15) << ErrorFrame << '\t'
                 << setw(20) << ErrorBits << '\t' << setw(20) << FER << '\t' << setw(20) << BER << '\t' << setw(15)
                 << LT3ErrBitFrame << '\t' << endl;

            //Output the noise seed to facilitate breakpoint continuation
            tout << "const unsigned long lastSeed[" << MAX_THREADS << "][3] = {\n";
            for (i = 0; i < MAX_THREADS; ++i) {
                tout << setw(4) << '{' << simulate[i].channel->RS.IX << "," << simulate[i].channel->RS.IY << ','
                     << simulate[i].channel->RS.IZ << "},\n";
            }
            tout << "}\n";
            tout.close();

            if (TestFrame > 1000 && ErrorFrame > 20) {
                break;
            }
            cout << setw(5) << snr << setw(20) << TestFrame << setw(15) << ErrorFrame << setw(20) << ErrorBits
                 << setw(20) << FER << setw(20) << BER << setw(15) << LT3ErrBitFrame << '\r' << flush;
        }

        gettimeofday(&t_end, NULL);
        double total_time = (t_end.tv_sec - t_start.tv_sec) + (t_end.tv_usec - t_start.tv_usec) / 1000000.0;
        cout << setw(5) << snr << setw(20) << TestFrame << setw(15) << ErrorFrame << setw(20) << ErrorBits << setw(20)
             << FER << setw(20) << BER << setw(15) << LT3ErrBitFrame << setw(15) << total_time << endl;
        fout << setw(5) << snr << '\t' << setw(20) << TestFrame << '\t' << setw(15) << ErrorFrame << '\t' << setw(20)
             << ErrorBits << '\t' << setw(20) << FER << '\t' << setw(20) << BER << '\t' << setw(15) << LT3ErrBitFrame
             << '\t' << setw(15) << total_time << '\t' << endl;
        fout.close();
        demodout << setw(5) << snr << '\t' << setw(20) << ModFER << '\t' << setw(20) << ModBER << '\t' << setw(20)
                 << ModSER << '\t' << endl;

        demodout.close();
    }
    getchar();
    return 0;
}
