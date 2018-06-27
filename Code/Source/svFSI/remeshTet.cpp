/* Copyright (c) Stanford University, The Regents of the University of
 *               California, and others.
 *
 * All Rights Reserved.
 *
 * See Copyright-SimVascular.txt for additional details.
 *
 * Permission is hereby granted, free of charge, to any person obtaining
 * a copy of this software and associated documentation files (the
 * "Software"), to deal in the Software without restriction, including
 * without limitation the rights to use, copy, modify, merge, publish,
 * distribute, sublicense, and/or sell copies of the Software, and to
 * permit persons to whom the Software is furnished to do so, subject
 * to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included
 * in all copies or substantial portions of the Software.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS
 * IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED
 * TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A
 * PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER
 * OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
 * EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
 * PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
 * PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
 * LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
 * NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
 * SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 */

    #include "tetgen.h"
    #include <iostream>
    #include <stdio.h>
    #include <stdlib.h>
    
    class tetOptions {
      public:
         double maxRadRatio;
         double minDihedAng;
         double maxEdgeSize;
         int optimLevel;
         int optimScheme;
         tetOptions();
         double maxTetVol(double r) {return (pow(r,3)/(6.0*sqrt(2)));}
   };
               
   tetOptions::tetOptions () {
      maxRadRatio = 1.15;
      minDihedAng = 10.0;
      maxEdgeSize = 0.15;
      optimLevel  = 2;
      optimScheme = 7;
   }

   extern"C"
   {
      void remesh3d_tetgen_(const int* nPoints, const int* nFacets, const double* pointList, \
                       const int* facetList, const double* params, int* pOK)
      {
         tetgenio in, out;
         tetgenio::facet *f;
         tetgenio::polygon *p;
         char fname [250];
         char switches [250];
         tetOptions options;
         
         *pOK = 0;
         in.firstnumber = 1;
         in.numberofpoints = *nPoints;
         in.pointlist = new REAL [in.numberofpoints * 3];
         for (int i=0; i < in.numberofpoints; i++)
         {
            for (int j=0; j < 3; j++) {
               in.pointlist[3*i+j] = *pointList;
               ++pointList;
            }
         }
         
         in.numberoffacets = *nFacets;
         in.facetlist = new tetgenio::facet[in.numberoffacets];
         in.facetmarkerlist = new int[in.numberoffacets];
         
         for (int i=0; i < in.numberoffacets; i++)
         {
            f = &in.facetlist[i];
            f->numberofpolygons = 1;
            f->polygonlist = new tetgenio::polygon[f->numberofpolygons];
            f->numberofholes = 0;
            f->holelist = NULL;
            p = &f->polygonlist[0];
            p->numberofvertices = 3;
            p->vertexlist = new int [p->numberofvertices];
            for (int j=0; j < 3; j++)
            {
               p->vertexlist[j] = *facetList;
               ++facetList;
            }
            in.facetmarkerlist[i] = 0;
         }
         
         options.maxRadRatio = *params;
         options.minDihedAng = *(++params);
         options.maxEdgeSize = *(++params);
         
         std::cout << " Using parameter <maxRadRatio> " << options.maxRadRatio << "\n";
         std::cout << " Using parameter <minDihedAng> " << options.minDihedAng << "\n";
         std::cout << " Using parameter <maxEdgeSize> " << options.maxEdgeSize << "\n\n";

         int len;
         len = sprintf(switches,"pYq%.2f/%.1fa%8.3eO%d/%d",    \
               options.maxRadRatio, options.minDihedAng,   \
               options.maxTetVol(options.maxEdgeSize),      \
               options.optimLevel,   options.optimScheme);
         
         if ( len > 250 )
         {
            std::cout << "    ERROR: Length of switch exceeded limit (250 char)\n";
            *pOK = -1;
            return;
         }
         tetrahedralize(switches, &in, &out);
         
         strcpy(fname, "new-vol-mesh");
         out.save_nodes(fname);
         out.save_elements(fname);
         
         return;
      }
   }
