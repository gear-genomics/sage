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

}

#endif
