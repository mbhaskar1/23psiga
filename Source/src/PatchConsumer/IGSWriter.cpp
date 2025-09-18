/* copyright(c)Jorg Peters [jorg.peters@gmail.com] */

#include "IGSWriter.hpp"

IGSWriter::IGSWriter(const std::string a_OutFile)
{
    m_OutFile = fopen(a_OutFile.c_str(),"w");
}


/* write knots of bb-form in IGES */
int IGSWriter::knots(int a_Deg1, int a_Bbase, int a_Ffctr, FILE* fp, int a_PerLine)
{
    int	i,j;

    for (i=0; i< a_Deg1; i++)
    {
    	if ((i%a_PerLine==0) && (i != 0))
        {
    	    fprintf(fp,"%8dP%7d\n", a_Ffctr, a_Bbase);
    	    a_Bbase++;
    	}
    	fprintf(fp,"0.00000,");
    }
    for (; i< 2*a_Deg1; i++)
    {
    	if (i%a_PerLine==0)
        {
    	    fprintf(fp, "%8dP%7d\n", a_Ffctr, a_Bbase);
    	    a_Bbase++;
	    }
	    fprintf(fp,"1.00000,");
    }
    j = (2*a_Deg1) % a_PerLine;
    if (j != 0)
    {
	    j = a_PerLine - j;
	    for (i=0; i < j; i++)
        {
		  fprintf(fp,"        ");
	    }
    }

    fprintf(fp, "%8dP%7d\n", a_Ffctr,a_Bbase);

    a_Bbase++;

    return(a_Bbase);
}


/*
 *  Destructor.
 *  All it does is close the file.
 */
IGSWriter::~IGSWriter()
{
    fprintf(m_OutFile, "copyright(c)Jorg Peters [jorg.peters@gmail.com]                         S      1\n");

    int k, k1, w, w1, rows, flen[4], i, m, j, col, cols, fc, sds, dg;
    double	h;

    //int facelen = 2;
    int t_Bbctr = 1;
    int t_Ffctr = 1;

    for (fc=0; fc < m_Patches.size(); fc++)
    {

        auto t_Patch = m_Patches[fc];
        int t_DegU = t_Patch.m_DegU;
        int t_DegV = t_Patch.m_DegV;

        //k = (dg+1)*2;
        k = (t_DegU + 1) + (t_DegV + 1);
        k1 = (k % 8 != 0) + k / 8;  /* knots */
        w = (t_DegU + 1) * (t_DegV + 1);
        w1 = (w % 8 != 0) + w / 8;  /* weights */

        /* deg-line, knots, weights, xyz location, param */
        rows  = (1 + 2 * k1 + w1 + w + 1);  /* size of one block */

        fprintf(m_OutFile,"     128%8d       0       1       0       0       0        00000000D%7d\n", t_Bbctr, t_Ffctr);
        fprintf(m_OutFile,"     128%8d       8      %d       2                NurbSurf       0D%7d\n", 0, rows, t_Ffctr+1);

        t_Bbctr += rows;
        t_Ffctr += 2;
    }
    flen[0] = 1;
    flen[1] = 10;

    /* data entries */
    t_Ffctr = 1; t_Bbctr = 1;
    for (fc=0; fc < m_Patches.size(); fc++)
    {

        auto t_Patch = m_Patches[fc];
        int t_DegU = t_Patch.m_DegU;
        int t_DegV = t_Patch.m_DegV;

        /* HEADER */
        fprintf(m_OutFile, "128,%7d,%7d,%7d,%7d,0,0,1,0,0,%26dP%7d\n",
                t_DegV, t_DegU, t_DegV, t_DegU, t_Ffctr, t_Bbctr);
        t_Bbctr++;

        /* KNOTS */
        w = (t_DegU + 1) * (t_DegV + 1);
        t_Bbctr = IGSWriter::knots(t_DegV + 1, t_Bbctr, t_Ffctr, m_OutFile, 8);
        t_Bbctr = IGSWriter::knots(t_DegU + 1, t_Bbctr, t_Ffctr, m_OutFile, 8);

        /* the w==RATIONAL coordinate */
        cols = 8;
        col = 0;
        for (i=0; i<= t_DegU; i++)
        {
            for (j=0; j<= t_DegV; j++)
            {
                if ((col % cols==0) && (col != 0)) {  /* typeset */
                    fprintf(m_OutFile, "%8dP%7d\n", t_Ffctr,t_Bbctr);
                    t_Bbctr++;
                }
                h = 1;
                fprintf(m_OutFile, "%7.5f,", h);
                col++;
            }
        }

        /* finish line */
        col = col%cols;
        if (col != 0)
        {
            col = cols-col;
            for (i=0; i< col; i++)
            {
                fprintf(m_OutFile,"        ");
            }
        }

        fprintf(m_OutFile, "%8dP%7d\n", t_Ffctr, t_Bbctr);
        t_Bbctr++;

        /* the XYZ coordinates */
        for (i=0; i<= t_DegU; i++)
        {
            for (j=0; j<= t_DegV; j++)
            {
                fprintf(m_OutFile,"%20e,%20e,%20e,%9dP%7d\n",
                        t_Patch.m_BBcoefs[0](i, j), t_Patch.m_BBcoefs[1](i, j), t_Patch.m_BBcoefs[2](i, j), t_Ffctr,t_Bbctr);
                t_Bbctr++;
            }
        }

        fprintf(m_OutFile, "0.00000,1.00000,0.00000,1.00000;%40dP%7d\n", t_Ffctr, t_Bbctr);
        t_Bbctr++;
        t_Ffctr += 2;
    }
    flen[2] = t_Ffctr;
    flen[3] = t_Bbctr;

    /* structure of file */
    fprintf(m_OutFile, "S%7dG%7dD%7dP%7d%40dT%7d\n", flen[0], flen[1], flen[2] - 1, flen[3] - 1, 1, 1);

    fclose(m_OutFile);
}

/*
 *  Calls WriteFromQueue in a new thread
 *  and then detaches that thread.
 */
void IGSWriter::start()
{
    return;
}

/*
 *  Signals to WriteFromQueue to stop once
 *  it has processed all the patches in
 *  the queue.
 */
void IGSWriter::stop()
{
    return;
}

/*
 *  Adds a patch to the queue in a
 *  thread-safe way.
 */
void IGSWriter::consume(const Patch a_Patch)
{
    writePatch(a_Patch);
    return;
}

void IGSWriter::writePatch(const Patch a_Patch)
{
    // Check that the patch has the correct
    // number of control points for its degree
    if(!a_Patch.isValid())
    {
        std::cout << "Patch is not valid,  write is aborting!\n";
        return;
    }

    // Save the patch into buffer
    m_Patches.push_back(a_Patch);

    return;
}
