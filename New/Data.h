#ifndef M31Globulars_Data_h
#define M31Globulars_Data_h

#include <vector>

class Data
{
    private:
        std::vector<double> x, y, v, verr, metallicity;        

    public:
        Data(const char* filename);

};

#endif
