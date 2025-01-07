#include <sstream>
#include <string>
#include <limits>
#include <fstream>
#include <vector>
#include <chrono>
#include <iomanip>
#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include <pybind11/stl.h>
#include <vector>
#include <algorithm>

// https://pybind11.readthedocs.io/en/latest/classes.html#instance-and-static-fields

using namespace std;
namespace py = pybind11;

struct GENEAData {
    std::vector<float> accelerationX;
    std::vector<float> accelerationY;
    std::vector<float> accelerationZ;
    float samplingRate;
    std::string measurementStartTime;

    const py::array getAccelerationX() const { return py::cast(accelerationX); }
    const py::array getAccelerationY() const { return py::cast(accelerationY); }
    const py::array getAccelerationZ() const { return py::cast(accelerationZ); }
    const float getSamplingRate() const { return samplingRate; }
    const py::str getMeasurementStartTime() const { return py::cast(measurementStartTime); }
};

static std::string pageTimeToStr(string& str) {
    str = str.replace(str.rfind(':'), 1, ".").substr(str.find(":") + 1);
    return str;
}

int parseBinFileHeader(std::istream& input_file, int fileHeaderSize, int linesToAxesCalibration,
    double(&gainVals)[3], int(&offsetVals)[3]) {
    // skip the first i lines in the file
    auto max_streamsize = std::numeric_limits<std::streamsize>::max();
    for (int i = 0; i < linesToAxesCalibration; i++) {
        input_file.ignore(max_streamsize, '\n');
    }
    // read axes calibration lines for gain and offset values
    // data like -> x gain:25548 \n x offset:574 ... Volts:300 \n Lux:800
    std::string line;
    input_file.ignore(max_streamsize, ':');
    input_file >> gainVals[0]; // xGain
    input_file.ignore(max_streamsize, ':');
    input_file >> offsetVals[0]; // xOffset
    input_file.ignore(max_streamsize, ':');
    input_file >> gainVals[1]; // y
    input_file.ignore(max_streamsize, ':');
    input_file >> offsetVals[1]; // y
    input_file.ignore(max_streamsize, ':');
    input_file >> gainVals[2]; // z
    input_file.ignore(max_streamsize, ':');
    input_file >> offsetVals[2]; // z

    int volts, lux, numBlocksTotal;

    input_file.ignore(max_streamsize, ':');
    input_file >> volts; // volts
    input_file.ignore(max_streamsize, ':');
    input_file >> lux; // lux

    input_file.ignore(max_streamsize, '\n'); // 9 blank
    input_file.ignore(max_streamsize, '\n'); // 10 memory status header

    input_file.ignore(max_streamsize, ':');
    input_file >> numBlocksTotal; // 11
    input_file.ignore(max_streamsize, '\n');

    // ignore remaining header lines in bin file
    for (int i = 0; i < fileHeaderSize - linesToAxesCalibration - 11; i++) {
        input_file.ignore(max_streamsize, '\n');
    }
    return numBlocksTotal;
}

int getSignedIntFromHex(const std::string& hex) {
    // input hex base is 16
    int rawVal = std::stoll(hex, nullptr, 16);
    int unsignedLimit = 4096; // 2^[length*4] #i.e. 3 hexBytes (12 bits)
    // limit = 4096
    int signedLimit = 2048; // 2^[length*(4-1)] #i.e. 3 hexBytes - 1 bit (11
    // bits) limit = 2048
    if (rawVal >= signedLimit) {
        rawVal = rawVal - unsignedLimit;
    }
    return rawVal;
}

int getLight(const std::string& hex) {
    // input hex base is 16
    int rawVal = std::stoll(hex, nullptr, 16);
    return rawVal;
}

// From https://github.com/wadpac/GGIRread
GENEAData GENEActivReader(std::string filename) {
    std::size_t start = 0;
    std::size_t end = 0;
    bool progress_bar = false;
    int fileHeaderSize = 59;
    int linesToAxesCalibration = 47;
    int blockHeaderSize = 9;
    int statusOK = -1;
    double sampleRate = -1;
    int errCounter = 0;

    std::vector<long> time_array;
    std::vector<float> x_array, y_array, z_array, temperature_array, lux_array;

    auto max_streamsize = std::numeric_limits<std::streamsize>::max();

    std::size_t numBlocksTotal = 0;
    double freq = 0.0;
    std::string measurementStartTime = "0";

    try {
        std::ifstream input_file(filename);
        // Read header to determine mfrGain and mfrOffset values
        double mfrGain[3];
        int mfrOffset[3];
        int numBlocksTotalint = parseBinFileHeader(input_file, fileHeaderSize, linesToAxesCalibration, mfrGain, mfrOffset);
        numBlocksTotal = numBlocksTotalint;

        std::size_t blockCount = 0;
        std::string header;
        long blockTime = 0;  // Unix millis
        long lastvalue = 0;  // Unix millis
        double temperature = 0.0;

        std::string data;
        std::string timeFmtStr = "Page Time:%Y-%m-%d %H:%M:%S:";

        if (start == 0) {
            start = 1;
        }
        if (end == 0) {
            end = numBlocksTotal;
        }
        std::string line;
        while (std::getline(input_file, line)) {
            ++blockCount;
            if (blockCount >= start && blockCount <= end) {
                // header: "Recorded Data" (0), serialCode (1), seq num (2),
                // blockTime (3), unassigned (4), temp (5), batteryVolt (6),
                // deviceStatus (7), freq (8), data (9)
                for (int i = 1; i < blockHeaderSize; i++) {
                    try {
                        std::getline(input_file, header);
                        if (i == 3) {
                            std::stringstream ss(header);
                            int milliseconds;
                            ss >> milliseconds;
                            blockTime = lastvalue + milliseconds;
                            if (blockCount == 1) {
                                measurementStartTime = pageTimeToStr(header);
                            }
                        }
                        else if (i == 5) {
                            std::stringstream ss(header);
                            ss.ignore(max_streamsize, ':');
                            ss >> temperature;
                        }
                        else if (i == 8) {
                            std::stringstream ss(header);
                            ss.ignore(max_streamsize, ':');
                            ss >> freq;
                        }
                    }
                    catch (const std::exception& e) {
                        errCounter++;
                        continue;
                    }
                }
                sampleRate = freq;

                // now process hex data
                std::getline(input_file, data);

                // raw reading values
                std::size_t hexPosition = 0;
                int xRaw = 0;
                int yRaw = 0;
                int zRaw = 0;
                int lux = 0;
                double x = 0.0;
                double y = 0.0;
                double z = 0.0;
                double t = 0.0;

                int i = 0;
                while (hexPosition < data.size() - 1) {
                    try {
                        std::stringstream ss;
                        xRaw = getSignedIntFromHex(data.substr(hexPosition, 3));
                        yRaw = getSignedIntFromHex(data.substr(hexPosition + 3, 3));
                        zRaw = getSignedIntFromHex(data.substr(hexPosition + 6, 3));
                        // get last 10th and 11th heximal places, and convert to binary
                        lux = getLight(data.substr(hexPosition + 9, 3)) / 4;
                        // Update values to calibrated measure (taken from GENEActiv manual)
                        x = (xRaw * 100. - mfrOffset[0]) / mfrGain[0];
                        y = (yRaw * 100. - mfrOffset[1]) / mfrGain[1];
                        z = (zRaw * 100. - mfrOffset[2]) / mfrGain[2];

                        t = (double)blockTime + (double)i * (1.0 / freq) * 1000;  // Unix millis
                        lastvalue = t;
                        time_array.push_back(t);
                        x_array.push_back(x);
                        y_array.push_back(y);
                        z_array.push_back(z);
                        temperature_array.push_back(temperature);
                        lux_array.push_back(lux);
                        hexPosition += 12;
                        i++;
                    }
                    catch (const std::exception& ex) {
                        errCounter++;
                        std::string err_text = ex.what();
                        break;  // rest of this block could be corrupted
                    }
                }
            }
            else if (blockCount < start) {
                // skip this block
                for (int i = 1; i < blockHeaderSize; i++) {  // header
                    input_file.ignore(max_streamsize, '\n');
                }
                input_file.ignore(max_streamsize, '\n');     // hexdata
            }
            else {
                // after end, no need to scan further
                break;
            }

        }
        statusOK = 1;
    }
    catch (const std::exception& e) {
        statusOK = 0;
    }

    GENEAData data;
    data.accelerationX = x_array;
    data.accelerationY = y_array;
    data.accelerationZ = z_array;
    data.measurementStartTime = measurementStartTime;
    data.samplingRate = freq;
    return data;
}


PYBIND11_MODULE(pybindRetTest, m) {
    py::class_<GENEAData>(m, "GENEAData")
        .def("getAccelerationX", &GENEAData::getAccelerationX)
        .def("getAccelerationY", &GENEAData::getAccelerationY)
        .def("getAccelerationZ", &GENEAData::getAccelerationZ)
        .def("getSamplingRate", &GENEAData::getSamplingRate)
        .def("getMeasurementStartTime", &GENEAData::getMeasurementStartTime);

    m.def("GENEActivReader", &GENEActivReader, R"pbdoc(
        Reads GENEActiv data.
    )pbdoc");

#ifdef VERSION_INFO
    m.attr("__version__") = VERSION_INFO;
#else
    m.attr("__version__") = "dev";
#endif
}
