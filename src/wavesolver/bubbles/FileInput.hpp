#ifndef _FILE_INPUT_HPP
#define _FILE_INPUT_HPP

#include <cmath>
#include <vector>
#include <string>

#include <boost/algorithm/string.hpp>

typedef struct
{
    int bubNum;
    double freq; // Negative for bad/dead bubbles
    std::complex<double> xfrIn; // TODO: unnecessary here
} BubbleInputInfo;

static std::vector<BubbleInputInfo>
parseFreqFile(const std::string &fName)
{
    using namespace std;
    vector<BubbleInputInfo> bubInfo;

    ifstream in(fName.c_str());

    string line;

    bubInfo.clear();

    // First line is time
    getline(in, line);

    // Read first real line
    getline(in, line);

    vector<string> data;
    while (!line.empty())
    {
        boost::split(data, line, boost::is_any_of(" \n"), boost::token_compress_on);

        if (data.size() < 2)
        {
            break;
        }

        BubbleInputInfo curInfo;
        curInfo.bubNum = atoi(data[0].c_str());

        if (data.size() == 2)
        {
            // Bad bubble
            curInfo.freq = -1;
        }
        else
        {
            // Good bubble

            curInfo.freq = atof(data[1].c_str());

            vector<string> xfrData;
            boost::split(xfrData, data[2], boost::is_any_of("(),"), boost::token_compress_on);

            curInfo.xfrIn = complex<double>(atof(xfrData[1].c_str()), atof(xfrData[2].c_str()));
        }

        bubInfo.push_back(curInfo);

        getline(in, line);
    }

    return bubInfo;
}

#endif // _FILE_INPUT_HPP

