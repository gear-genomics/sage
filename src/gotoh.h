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

#ifndef GOTOH_H
#define GOTOH_H

#include <iostream>
#include "align.h"

namespace sage
{
  
  template<typename TProfile, typename TAlign, typename TAlignConfig, typename TScoreObject>
  inline int
  gotoh(TProfile const& p1, TProfile const& p2, TAlign& align, TAlignConfig const& ac, TScoreObject const& sc) {
    typedef typename TScoreObject::TValue TScoreValue;

    // DP Matrix
    typedef boost::multi_array<TScoreValue, 2> TMatrix;
    std::size_t m = p1.shape()[1];
    std::size_t n = p2.shape()[1];
    TMatrix s(boost::extents[m+1][n+1]);
    TMatrix h(boost::extents[m+1][n+1]);
    TMatrix v(boost::extents[m+1][n+1]);

    // Initialization
    for(std::size_t col = 1; col <= n; ++col) {
      v[0][col] = -sc.inf;
      s[0][col] = _horizontalGap(ac, 0, m, sc.go + col * sc.ge);
      h[0][col] = _horizontalGap(ac, 0, m, sc.go + col * sc.ge);
    }
    for(std::size_t row = 1; row <= m; ++row) {
      h[row][0] = -sc.inf;
      s[row][0] = _verticalGap(ac, 0, n, sc.go + row * sc.ge);
      v[row][0] = _verticalGap(ac, 0, n, sc.go + row * sc.ge);
    }
    s[0][0] = 0;
    v[0][0] = -sc.inf;
    h[0][0] = -sc.inf;

    // Recursion
    for(std::size_t row = 1; row <= m; ++row) {
      for(std::size_t col = 1; col <= n; ++col) {
	h[row][col] = std::max(s[row][col-1] + _horizontalGap(ac, row, m, sc.go + sc.ge), h[row][col-1] + _horizontalGap(ac, row, m, sc.ge));
	v[row][col] = std::max(s[row-1][col] + _verticalGap(ac, col, n, sc.go + sc.ge), v[row-1][col] + _verticalGap(ac, col, n, sc.ge));
	s[row][col] = std::max(std::max(s[row-1][col-1] + _score(p1, p2, row-1, col-1, sc), h[row][col]), v[row][col]);
      }
    }

    // Trace-back
    std::size_t row = m;
    std::size_t col = n;
    char lastMatrix = 's';
    typedef std::vector<char> TTrace;
    TTrace trace;
    while ((row>0) || (col>0)) {
      if (lastMatrix == 's') {
	if (s[row][col] == h[row][col]) lastMatrix = 'h';
	else if (s[row][col] == v[row][col]) lastMatrix = 'v';
	else {
	  --row;
	  --col;
	  trace.push_back('s');
	}
      } else if (lastMatrix == 'h') {
	if (h[row][col] != h[row][col-1] + _horizontalGap(ac, row, m, sc.ge)) lastMatrix = 's';
	--col;
	trace.push_back('h');
      } else if (lastMatrix == 'v') {
	if (v[row][col] != v[row-1][col] + _verticalGap(ac, col, n, sc.ge)) lastMatrix = 's';
	--row;
	trace.push_back('v');
      }
    }

    // Create alignment
    _createAlignment(trace, p1, p2, align);

    // Score
    return s[m][n];
  }

}

#endif
