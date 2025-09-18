/* copyright(c)Jorg Peters [jorg.peters@gmail.com] */

#pragma once

#include <string>
#include <mutex>
#include <fstream>
#include <queue>

#include "PatchConsumer.hpp"
#include "../Patch/Patch.hpp"

class BVWriter : public PatchConsumer
{
public:
    BVWriter(const std::string a_OutFile);
    ~BVWriter();
    void start();
    void stop();
    void consume(const Patch a_Patch);
private:
    std::ofstream m_OutFile;
    void writePatch(const Patch a_Patch);
};
