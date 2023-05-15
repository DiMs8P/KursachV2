#pragma once
#include "Config.h"
#include "vector"
using namespace std;

class TimeParser
{
public:
    TimeParser(TimeInfo1 timeInfo)
    {
        _timeInfo = timeInfo;
    }

    void Parse(vector<double>& output)
    {
        output.push_back(_timeInfo.StartTime);
        for (int i = 1; i < _timeInfo.TimesNum; i++)
        {
            output.push_back(output[i-1] + _timeInfo.InitialStep * pow(_timeInfo.StepMultiplier, i - 1));
        }
    }
private:
    TimeInfo1 _timeInfo;
};
