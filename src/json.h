/*
============================================================================
Sage: Sanger Trace Alignment
============================================================================
Copyright (C) 2017 Tobias Rausch

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.
============================================================================
Contact: Tobias Rausch (rausch@embl.de)
============================================================================
*/

#ifndef JSON_H
#define JSON_H

#include <boost/progress.hpp>
#include "abif.h"

namespace sage
{

  template<typename TStream, typename TAlign>
  inline void
  _traceAlignJsonOut(TStream& rfile, ReferenceSlice const& rs, TAlign const& align) {
    rfile << "\"refchr\": \"" << rs.chr << "\"," << std::endl;
    rfile << "\"refpos\": " << rs.pos << "," << std::endl;
    rfile << "\"refalign\": \"";
    for(uint32_t j = 0; j<align.shape()[1]; ++j) rfile << align[1][j];
    rfile << "\"," << std::endl;
    rfile << "\"queryalign\": \"";
    for(uint32_t j = 0; j<align.shape()[1]; ++j) rfile << align[0][j];
    rfile << "\"," << std::endl;
    rfile << "\"forward\": " << rs.forward << ',' << std::endl;
  }

  template<typename TAlign>
  inline void
  traceAlignJsonOut(std::string const& outfile, teal::BaseCalls& bc, teal::Trace const& tr, ReferenceSlice const& rs, TAlign const& align) {
    // Output trace
    std::ofstream rfile(outfile.c_str());
    rfile << "{" << std::endl;
    _traceAlignJsonOut(rfile, rs, align);
    teal::_traceJsonOut(rfile, bc, tr);
    rfile << "}" << std::endl;
    rfile.close();  
  }

  template<typename TAlign>
  inline void
  gappyFuct( teal::BaseCalls& nbc,  teal::Trace& ntr, teal::BaseCalls& bc, teal::Trace const& tr, TAlign const& align) {
    int abioffset = 49;
   
    typedef teal::Trace::TMountains TMountains; 
    typedef teal::Trace::TValue TValue;
    typedef std::vector<uint32_t> TVecVal;
    TVecVal insPos;
    TVecVal insSize; 

    // Calculate average basecall in this file
    uint32_t step = 6;
    if (bc.bcPos.size() > 151) {
        step = (uint32_t) ((float) bc.bcPos[150] - (float) bc.bcPos[50] ) / 200.0 + 0.9;
    }
    // Find the gap positions 
    uint32_t pos = abioffset;
    char mod = 's';
    uint32_t gapsize = 0;    
    for(uint32_t j = 0; j<align.shape()[1]; ++j) {
        if (align[0][j] == '-') {
            if (mod != 's') {
                gapsize++;
            }
            mod = '-';
        } else {
            if (mod == '-') {
                if (pos + 1 < bc.bcPos.size()) {
                    uint32_t insert = (uint32_t) (( bc.bcPos[pos + 1] + bc.bcPos[pos] ) / 2.0);
                    insPos.push_back(insert);
                    insSize.push_back(gapsize);
                }
            }
            gapsize = 0;
            mod = 'b';
            pos++;
        }
    }
    // Rewrite the arrays
    uint32_t offset = 0;
    int curPos = -1;
    if (curPos + 1 < insPos.size()) {
        curPos++;
    }
    uint32_t bcpos = 0;
    TValue idx = bc.bcPos[0];
    ntr.traceACGT.push_back(TMountains());   
    ntr.traceACGT.push_back(TMountains());
    ntr.traceACGT.push_back(TMountains());
    ntr.traceACGT.push_back(TMountains());
 
    for(uint32_t i = 0; i<tr.traceACGT[0].size(); ++i) {
        ntr.traceACGT[0].push_back(tr.traceACGT[0][i]);        
        ntr.traceACGT[1].push_back(tr.traceACGT[1][i]);
        ntr.traceACGT[2].push_back(tr.traceACGT[2][i]);
        ntr.traceACGT[3].push_back(tr.traceACGT[3][i]);
        if (idx == i) {
            nbc.bcPos.push_back(idx + offset);
            nbc.primary.push_back(bc.primary[bcpos]);
            nbc.secondary.push_back(bc.secondary[bcpos]);
            if (bcpos < bc.bcPos.size() - 1) idx = bc.bcPos[++bcpos];
        }
        if ((curPos > -1) && (i == insPos[curPos])) {
            for(uint32_t k = 0; k < insSize[curPos]; ++k) {
                for(uint32_t n = 0; n < step; ++n) {
                    offset++;
                    ntr.traceACGT[0].push_back(-99);
                    ntr.traceACGT[1].push_back(-99);
                    ntr.traceACGT[2].push_back(-99);
                    ntr.traceACGT[3].push_back(-99);
                }
                nbc.bcPos.push_back(i + offset);
                nbc.primary.push_back('-');
                nbc.secondary.push_back('-');
                for(uint32_t n = 0; n < step; ++n) {
                    offset++;
                    ntr.traceACGT[0].push_back(-99);
                    ntr.traceACGT[1].push_back(-99);
                    ntr.traceACGT[2].push_back(-99);
                    ntr.traceACGT[3].push_back(-99);
                }
            }
            ntr.traceACGT[0].push_back(tr.traceACGT[0][i]);
            ntr.traceACGT[1].push_back(tr.traceACGT[1][i]);
            ntr.traceACGT[2].push_back(tr.traceACGT[2][i]);
            ntr.traceACGT[3].push_back(tr.traceACGT[3][i]);
            if (curPos + 1 < insPos.size()) {
                curPos++;
            }
        }
    }
  }


}

#endif
