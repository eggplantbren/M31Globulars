#ifndef M31Globulars_Data_h
#define M31Globulars_Data_h

#include <vector>

class Data
{
    private:
        std::vector<double> x, y, v, verr, metallicity;  
        std::vector<double> theta;      

    public:
        Data(const char* filename);

        friend class MyModel;
};

#endif
