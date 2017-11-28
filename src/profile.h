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

#ifndef PROFILE_H
#define PROFILE_H

#include <boost/progress.hpp>

namespace sage
{


  template<typename TProfile>
  inline void
  createProfile(teal::Trace const& tr, teal::BaseCalls const& bc, TProfile& p) {
    p.resize(boost::extents[6][bc.bcPos.size()]);   // 'A', 'C', 'G', 'T', 'N', '-'
    for(uint32_t j = 0; j < bc.bcPos.size(); ++j) {
      double totalsig = 0;
      for(uint32_t k = 0; k<4; ++k) totalsig += tr.traceACGT[k][bc.bcPos[j]];
      p[4][j] = 0;
      p[5][j] = 0;
      for(uint32_t k = 0; k<4; ++k) p[k][j] = ((double) (tr.traceACGT[k][bc.bcPos[j]]) / totalsig);
    }
  }
  
  template<typename TProfile>
  inline void
  createProfile(teal::Trace const& tr, teal::BaseCalls const& bc, TProfile& p, int32_t const trimleft, int32_t const trimright) {
    if (trimleft + trimright >= (int32_t) bc.bcPos.size()) return createProfile(tr, bc, p);
    else {
      int32_t sz = bc.bcPos.size() - (trimleft + trimright);
      p.resize(boost::extents[6][sz]);   // 'A', 'C', 'G', 'T', 'N', '-'
      for(int32_t j = trimleft; j < (trimleft + sz); ++j) {
	double totalsig = 0;
	for(uint32_t k = 0; k<4; ++k) totalsig += tr.traceACGT[k][bc.bcPos[j]];
	p[4][j-trimleft] = 0;
	p[5][j-trimleft] = 0;
	for(uint32_t k = 0; k<4; ++k) p[k][j-trimleft] = ((double) (tr.traceACGT[k][bc.bcPos[j]]) / totalsig);
      }
    }
  }

  template<typename TProfile>
  inline void
  createProfile(std::string const& s, TProfile& p) {
    typedef typename TProfile::index TPIndex;
    p.resize(boost::extents[6][s.size()]);   // 'A', 'C', 'G', 'T', 'N', '-'
    for (std::size_t j = 0; j < s.size(); ++j) {
      for(TPIndex k = 0; k < 6; ++k) p[k][j] = 0;
      if ((s[j] == 'A') || (s[j] == 'a')) p[0][j] += 1;
      else if ((s[j] == 'C') || (s[j] == 'c')) p[1][j] += 1;
      else if ((s[j] == 'G') || (s[j] == 'g')) p[2][j] += 1;
      else if ((s[j] == 'T') || (s[j] == 't')) p[3][j] += 1;
      else if ((s[j] == 'N') || (s[j] == 'n')) p[4][j] += 1;
      else if (s[j] == '-') p[5][j] += 1;
    }
  }

  template<typename TChar, typename TProfile>
  inline void
  createProfile(boost::multi_array<TChar, 2> const& a, TProfile& p) {
    typedef typename boost::multi_array<TChar, 2>::index TAIndex;
    typedef typename TProfile::index TPIndex;
    p.resize(boost::extents[6][a.shape()[1]]);   // 'A', 'C', 'G', 'T', 'N', '-'
    for (TAIndex j = 0; j < (TAIndex) a.shape()[1]; ++j) {
      for(TPIndex k = 0; k < 6; ++k) p[k][j] = 0;
      int sum = 0;
      for(TAIndex i = 0; i < (TAIndex) a.shape()[0]; ++i) {
	++sum;
	if ((a[i][j] == 'A') || (a[i][j] == 'a')) p[0][j] += 1;
	else if ((a[i][j] == 'C') || (a[i][j] == 'c')) p[1][j] += 1;
	else if ((a[i][j] == 'G') || (a[i][j] == 'g')) p[2][j] += 1;
	else if ((a[i][j] == 'T') || (a[i][j] == 't')) p[3][j] += 1;
	else if ((a[i][j] == 'N') || (a[i][j] == 'n')) p[4][j] += 1;
	else if (a[i][j] == '-') p[5][j] += 1;
	else --sum;
      }
      for(TPIndex k = 0; k<6; ++k) p[k][j] /= sum;
    }
  }

}

#endif
