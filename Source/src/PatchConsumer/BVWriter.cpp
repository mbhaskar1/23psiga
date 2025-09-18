/* copyright(c)Jorg Peters [jorg.peters@gmail.com] */

#include <assert.h>
#include <iomanip>
#include <thread>
#include <cmath>

#include "BVWriter.hpp"

/*
 *  Constructor.
 *  All it does is open a file.
 */
BVWriter::BVWriter(const std::string a_OutFile)
{
    m_OutFile.open(a_OutFile);
}

/*
 *  Destructor.
 *  All it does is close the file.
 */
BVWriter::~BVWriter()
{
    m_OutFile.close();
}

/*
 *  Calls WriteFromQueue in a new thread
 *  and then detaches that thread.
 */
void BVWriter::start()
{
    /*
     *  Do nothing
     */
    return;
}

/*
 *  Signals to WriteFromQueue to stop once
 *  it has processed all the patches in
 *  the queue.
 */
void BVWriter::stop()
{
    /*
     *  Do nothing
     */
    return;
}

/*
 *  Adds a patch to the queue in a
 *  thread-safe way.
 */
void BVWriter::consume(const Patch a_Patch)
{
    writePatch(a_Patch);
    return;
}

/*
 *  Writes out a patch into a bv file with the correct format.
 *
 *  See https://www.cise.ufl.edu/research/SurfLab/bview/
 *  for information on the file format.
 */
void BVWriter::writePatch(Patch a_Patch)
{
    // Check that the patch has the correct
    // number of control points for its degree
    if(!a_Patch.isValid())
    {
        std::cout << "Patch is not valid,  write is aborting!\n";
        return;
    }

    // Patch Group
    m_OutFile << a_Patch.m_Group << std::endl;

    // Patch Type
    m_OutFile << a_Patch.m_PatchType << std::endl;

    // Degree in u and v directions
    const int t_NumOfFloatingDigit = 15;
    m_OutFile << std::fixed << std::setprecision(t_NumOfFloatingDigit) <<
        std::to_string(a_Patch.m_DegU) << " " << std::to_string(a_Patch.m_DegV) << std::endl;

    // Control points
    double fd = pow(10,t_NumOfFloatingDigit);
    std::vector<double> t_BB = {0, 0, 0};
    for(int row=0; row<a_Patch.m_DegU+1; ++row)
    {
        for(int column=0; column<a_Patch.m_DegV+1; ++column)
        {
            for(int i=0; i<t_BB.size(); i++)
            {
                t_BB[i] = round(a_Patch.m_BBcoefs[i](row, column) * fd);
                t_BB[i] = (long int)(t_BB[i])%2==0 ? t_BB[i]/fd : (t_BB[i]+1)/fd;
            }
            m_OutFile << t_BB[0] << " " << t_BB[1] << " " << t_BB[2] << std::endl;
        }
    }

    return;
}
