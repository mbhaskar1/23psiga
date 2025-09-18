/* copyright(c)Jorg Peters [jorg.peters@gmail.com] */

#include "ReadCSV2Matrix.hpp"
#include <iostream>
#include <fstream>
#include <string>

Matrix read_csv_as_matrix(const std::string a_File, const int a_Rows, const int a_Cols)
{
    std::ifstream t_File(a_File);
    std::string t_Line;

    int t_Row = 0;
    int t_Col = 0;

    Matrix t_Matrix(a_Rows, a_Cols);

    if (t_File.is_open())
    {
        while (std::getline(t_File, t_Line))
        {
            char *t_Ptr = (char *) t_Line.c_str();
            int t_Len = t_Line.length();

            // Skip comments
            if(t_Ptr[0] == '#')
            {
                continue;
            }

            t_Col = 0;

            char *t_Start = t_Ptr;
            for (int i=0; i < t_Len; i++)
            {
                if (t_Ptr[i] == ',')
                {
                    t_Matrix(t_Row, t_Col++) = atof(t_Start);
                    t_Start = t_Ptr + i + 1;
                }
            }
            t_Matrix(t_Row, t_Col) = atof(t_Start);

            t_Row++;
        }
        t_File.close();
    }
    else
    {
        std::cout << "Error opening file" << std::endl;
    }
    return t_Matrix;
}
