/* copyright(c)Jorg Peters [jorg.peters@gmail.com] */

#pragma once

#include <string>
#include <mutex>
#include <fstream>
#include <queue>

#include "PatchConsumer.hpp"
#include "../Patch/Patch.hpp"

class IGSWriter : public PatchConsumer
{
public:
    IGSWriter(const std::string a_OutFile);
    ~IGSWriter();
    void start();
    void stop();
    void consume(const Patch a_Patch);
private:
    std::vector<Patch> m_Patches;
    FILE* m_OutFile;
    int m_NumOfFaces = 0;
    int knots(int dg1, int bbase, int ffctr, FILE* fp, int per_line);
    void writePatch(const Patch a_Patch);
};
