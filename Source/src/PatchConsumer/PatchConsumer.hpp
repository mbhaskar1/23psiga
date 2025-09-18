/* copyright(c)Jorg Peters [jorg.peters@gmail.com] */

#pragma once

#include "../Patch/Patch.hpp"

/*
 *  Inheriting from this class allows us to change
 *  what we do with the patches. e.g. customize output format.
 *  We use BVWriter by default and provide IGSWriter as another example.
 *
 *  Other Example (commented code snippet): 
 *     compute the derivative of a BB-patch in the u direction
 *     by scaling, by the degree in u, the difference of BB-coefficients  [i][j] and [i-1][j]
 *     in a double loop over all BB-coefficients (that have an [i-1] neighbor)
 *
 *  std::vector<std::vector<double>> get_du(const Patch& a_Patch)
 *  {
 *      std::vector<std::vector<double>> t_Du;
 *      for(int i=1; i<a_Patch.m_DegU+1; i++)   // note: start at 1 so that [i-1] is well-defined
 *      {
 *          std::vector<double> t_DuRow;
 *          for(int j=0; j<a_Patch.m_DegV+1; j++)
 *          {
 *              // scaling, by the degree in u, the difference of BB-coefficients  [i][j] and [i-1][j]
 *              t_DuRow.push_back(a_Patch.m_DegU * (a_Patch.m_BBcoefs[i][j] - a_Patch.m_BBcoefs[i-1][j]));
 *          }
 *          t_Du.push_back(t_DuRow);
 *      }
 *      return t_Du;
 *  }
 *
 */
class PatchConsumer
{
public:
    virtual ~PatchConsumer() = default;
    virtual void start() = 0;
    virtual void stop() = 0;
    virtual void consume(Patch a_Patch) = 0;
};
