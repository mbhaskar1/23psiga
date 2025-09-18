/* copyright(c)Jorg Peters [jorg.peters@gmail.com] */

#pragma once

#include <iostream>
#include "../Helper/Helper.hpp"

#include <OpenMesh/Core/Mesh/PolyMesh_ArrayKernelT.hh>
#include <Eigen/Dense>

using Eigen::MatrixXd;
using Eigen::VectorXd;

typedef MatrixXd Matrix;
typedef VectorXd Vector;
typedef OpenMesh::PolyMesh_ArrayKernelT<>::Point Point;

struct Patch
{
    Patch() {};
    Patch(int a_BiDeg) : m_DegU(a_BiDeg), m_DegV(a_BiDeg)
    {
        for(int i = 0; i < 3; i++){
            m_BBcoefs.push_back(Matrix(a_BiDeg+1, a_BiDeg+1));
        }
        initBBcoefs();
    };

    Patch(int a_BiDeg, std::string a_Group) : m_DegU(a_BiDeg), m_DegV(a_BiDeg), m_Group(a_Group)
    {
        for(int i = 0; i < 3; i++){
            m_BBcoefs.push_back(Matrix(a_BiDeg+1, a_BiDeg+1));
        }
        initBBcoefs();
    };

    Patch(int a_DegU, int a_DegV) : m_DegU(a_DegU), m_DegV(a_DegV)
    {
        for(int i = 0; i < 3; i++){
            m_BBcoefs.push_back(Matrix(a_DegU+1, a_DegV+1));
        }
        initBBcoefs();
    };

    Patch(int a_DegU, int a_DegV, std::string a_Group) : m_DegU(a_DegU), m_DegV(a_DegV), m_Group(a_Group)
    {
        for(int i = 0; i < 3; i++){
            m_BBcoefs.push_back(Matrix(a_DegU+1, a_DegV+1));
        }
        initBBcoefs();
    };

    void operator=(const Patch a_Other)
    {
        m_DegU = a_Other.m_DegU;
        m_DegV = a_Other.m_DegV;
        m_Group = a_Other.m_Group;
        m_BBcoefs = a_Other.m_BBcoefs;
    }

    bool isValid() const
    {
        int t_OrderU = m_DegU + 1;
        int t_OrderV = m_DegV + 1;
        int t_ExpectedCpts = t_OrderU * t_OrderV;

        for(int idx = 0; idx < 3; idx++) {
            int t_NumCpts = m_BBcoefs[idx].rows() * m_BBcoefs[idx].cols();
            if(t_NumCpts != t_ExpectedCpts)
            {
                return false;
            }
        }

        return true;
    }

    // Init each cpts as {0,0,0}
    void initBBcoefs()
    {
        for(int idx = 0; idx < 3; idx++){
            for(int i = 0; i < m_BBcoefs[idx].rows(); i++){
                for(int j = 0; j < m_BBcoefs[idx].cols(); j++){
                    m_BBcoefs[idx](i, j) = 0;
                }
            }
        }
    }

    void degRaise()
    {
        for(int i = 0; i < 3; i++)
        {
            m_BBcoefs[i] = Helper::degRaise(m_BBcoefs[i]);
        }
        m_DegU = m_BBcoefs[0].rows() - 1;
        m_DegV = m_BBcoefs[0].cols() - 1;
    }

    // Please find patch type in
    // https://www.cise.ufl.edu/research/SurfLab/bview/#file-format
    const std::string m_PatchType = "5";

    int m_DegU;
    int m_DegV;
    std::string m_Group = "Group 0 default";

    // m_BBcoefs should be a 3 x m_DegU+1 x m_DegV+1 array
    std::vector<Matrix> m_BBcoefs;
};

static Patch points_mat_to_patch(const Matrix& a_PointMat)
{
    int t_PointsPerCol = pow(a_PointMat.rows(), 0.5);
    Patch t_Patch(t_PointsPerCol-1);
    for(int i=0; i<a_PointMat.rows(); i++)
    {
        for(int idx = 0; idx < 3; idx++)
        {
            t_Patch.m_BBcoefs[idx](i/t_PointsPerCol, i%t_PointsPerCol) = a_PointMat(i,idx);
        }
    }
    return t_Patch;
}

static Patch points_mat_to_patch(const int a_PatchDegU, const int a_PatchDegV, const Matrix& a_PointMat)
{
    Patch t_Patch(a_PatchDegU, a_PatchDegV);
    for(int i=0; i<a_PatchDegU+1; i++)
    {
        for(int j=0; j<a_PatchDegV+1; j++)
        {
            int t_Index = i * (a_PatchDegV+1) + j;
            for(int idx = 0; idx < 3; idx++)
            {
                t_Patch.m_BBcoefs[idx](i, j) = a_PointMat(t_Index,idx);
            }
        }
    }
    return t_Patch;
}

static Patch points_mat_to_patch(const int a_PatchDegU, const int a_PatchDegV, const std::string a_Group, const Matrix& a_PointMat)
{
    Patch t_Patch(a_PatchDegU, a_PatchDegV, a_Group);
    for(int i=0; i<a_PatchDegU+1; i++)
    {
        for(int j=0; j<a_PatchDegV+1; j++)
        {
            int t_Index = i * (a_PatchDegV+1) + j;
            for(int idx = 0; idx < 3; idx++)
            {
                t_Patch.m_BBcoefs[idx](i, j) = a_PointMat(t_Index,idx);
            }
        }
    }
    return t_Patch;
}

static Patch points_mat_to_patch(const int a_PatchDegU, const int a_PatchDegV, const std::string a_Group, const Matrix& a_PointMat, const int a_StartIndex)
{
    Patch t_Patch(a_PatchDegU, a_PatchDegV, a_Group);
    for(int i=0; i<a_PatchDegU+1; i++)
    {
        for(int j=0; j<a_PatchDegV+1; j++)
        {
            int t_Index = a_StartIndex + i * (a_PatchDegV+1) + j;
            for(int idx = 0; idx < 3; idx++)
            {
                t_Patch.m_BBcoefs[idx](i, j) = a_PointMat(t_Index,idx);
            }
        }
    }
    return t_Patch;
}

static std::vector<Patch> points_mat_to_patches(const int a_NumOfPatch, const std::string a_Group, const Matrix& a_PointMat)
{
    std::vector<Patch> t_Patches;
    const int t_NumOfCCPtsPerPatch = a_PointMat.rows() / a_NumOfPatch;
    const int t_Deg = sqrt(t_NumOfCCPtsPerPatch) - 1;
    for(int i=0; i<a_NumOfPatch; i++ )
    {
        t_Patches.push_back(points_mat_to_patch(t_Deg, t_Deg, a_Group, a_PointMat, i*t_NumOfCCPtsPerPatch));
    }

    return t_Patches;
}

/*
 * For the patch whose height and width are not the same
 */
static std::vector<Patch> points_mat_to_patches(const int a_PatchDegU, const int a_PatchDegV, const std::string a_Group, const Matrix& a_PointMat)
{
    std::vector<Patch> t_Patches;
    const int t_NumOfCCPtsPerPatch = (a_PatchDegU + 1) * (a_PatchDegV + 1);
    const int t_NumOfPatch = a_PointMat.rows() / t_NumOfCCPtsPerPatch;

    for(int i=0; i<t_NumOfPatch; i++ )
    {
        t_Patches.push_back(points_mat_to_patch(a_PatchDegU, a_PatchDegV, a_Group, a_PointMat, i*t_NumOfCCPtsPerPatch));
    }

    return t_Patches;
}
