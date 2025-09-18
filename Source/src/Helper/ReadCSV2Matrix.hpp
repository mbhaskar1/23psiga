/* copyright(c)Jorg Peters [jorg.peters@gmail.com] */

#pragma once

#include <Eigen/Dense>

using Eigen::MatrixXd;

typedef MatrixXd Matrix;

Matrix read_csv_as_matrix(const std::string a_File, const int a_Rows, const int a_Cols);
