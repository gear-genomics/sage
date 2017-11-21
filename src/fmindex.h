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

#ifndef FMINDEX_H
#define FMINDEX_H

#include <boost/progress.hpp>

namespace sage
{

  struct ReferenceSlice {
    bool forward;
    uint32_t filetype;   //0: *fa.gz, 1: *.fa, 2: *.ab1
    uint32_t kmersupport;
    uint32_t pos;
    std::string chr;
    std::string refslice;
  };
  

  template<typename TConfig, typename TFMIdx>
  inline bool
  loadFMIdx(TConfig const& c, ReferenceSlice& rs, TFMIdx& fm_index) {
    
    // What kind of reference?
    std::ifstream ifile(c.genome.c_str(), std::ios::binary | std::ios::in);
    if (ifile.is_open()) {
      char fcode[4];
      ifile.seekg(0);
      ifile.read(fcode, 4);
      if (((uint8_t)fcode[0] == (uint8_t)0x1f) && ((uint8_t)fcode[1] == (uint8_t)0x8b)) {
	// Gzipped fasta
	rs.filetype = 0;
	boost::filesystem::path op = c.genome.parent_path() / c.genome.stem();
	boost::filesystem::path outfile(op.string() + ".dump");
	std::string index_file = op.string() + ".fm9";

	// Load FM index
	boost::posix_time::ptime now = boost::posix_time::second_clock::local_time();
	std::cout << '[' << boost::posix_time::to_simple_string(now) << "] " << "Load FM-Index" << std::endl;
	if (!load_from_file(fm_index, index_file)) {
	  std::cout << '[' << boost::posix_time::to_simple_string(now) << "] " << "Build FM-Index" << std::endl;
	  // Dump fasta
	  bool firstSeq = true;
	  std::ofstream tmpout(outfile.string().c_str());
	  std::ifstream file(c.genome.string().c_str(), std::ios_base::in | std::ios_base::binary);
	  boost::iostreams::filtering_streambuf<boost::iostreams::input> dataIn;
	  dataIn.push(boost::iostreams::gzip_decompressor());
	  dataIn.push(file);
	  std::istream instream(&dataIn);
	  std::string line;
	  while(std::getline(instream, line)) {
	    if (line.find(">") == 0) {
	      if (!firstSeq) tmpout << std::endl;
	      else firstSeq = false;
	    } else {
	      tmpout << boost::to_upper_copy(line);
	    }
	  }
	  tmpout << std::endl;
	  file.close();
	  tmpout.close();
	  
	  now = boost::posix_time::second_clock::local_time();
	  std::cout << '[' << boost::posix_time::to_simple_string(now) << "] " << "Create FM-Index" << std::endl;
      
	  // Build index
	  construct(fm_index, outfile.string().c_str(), 1);
	  store_to_file(fm_index, index_file);
	  boost::filesystem::remove(outfile);
	}
      } else if (std::string(fcode) == "ABIF") {
	rs.filetype = 2;
	
	boost::posix_time::ptime now = boost::posix_time::second_clock::local_time();
	std::cout << '[' << boost::posix_time::to_simple_string(now) << "] " << "Load ab1 wildtype" << std::endl;
	teal::Trace wt;
	if (!teal::readab(c.genome.string(), wt)) return false;
	teal::BaseCalls wtbc;
	teal::basecall(wt, wtbc, c.pratio);
	rs.chr = "wildtype";
	rs.refslice = wtbc.primary;
	construct_im(fm_index, rs.refslice.c_str(), 1);
      }
      else if (fcode[0] == '>') {
	rs.filetype = 1;
	
	// Single FASTA file
	rs.chr = "";
	std::string tmpfasta = "";
	boost::posix_time::ptime now = boost::posix_time::second_clock::local_time();
	std::cout << '[' << boost::posix_time::to_simple_string(now) << "] " << "Load FASTA reference" << std::endl;
	std::ifstream fafile(c.genome.string().c_str());
	if (fafile.good()) {
	  std::string line;
	  while(std::getline(fafile, line)) {
	    if (!line.empty()) {
	      if (line[0] == '>') {
		if (!rs.chr.empty()) {
		  std::cerr << "Only single-chromosome FASTA files are supported. If you have a multi-FASTA file please use bgzip and index the FASTA file with samtools faidx!" << std::endl;
		  return false;
		}
		rs.chr = line.substr(1);
	      } else {
		tmpfasta += boost::to_upper_copy(line);
	      }
	    }
	  }
	  fafile.close();
	}
	// Check FASTA
	rs.refslice = "";
	for(uint32_t k = 0; k < tmpfasta.size(); ++k)
	  if ((tmpfasta[k] == 'A') || (tmpfasta[k] == 'C') || (tmpfasta[k] == 'G') || (tmpfasta[k] == 'T') || (tmpfasta[k] == 'N')) rs.refslice += tmpfasta[k];
	if (rs.refslice.size() != tmpfasta.size()) {
	  std::cerr << "FASTA file contains nucleotides != [ACGTN]." << std::endl;
	  return false;
	}
	construct_im(fm_index, rs.refslice.c_str(), 1);
      } else {
	std::cerr << "Couldn't recognize reference file format!" << std::endl;
	return false;
      }
    }
    ifile.close();
    return true;
  }


  inline void
  reverseComplement(std::string& sequence) {
    std::string rev = boost::to_upper_copy(std::string(sequence.rbegin(), sequence.rend()));
    std::size_t i = 0;
    for(std::string::iterator revIt = rev.begin(); revIt != rev.end(); ++revIt, ++i) {
      switch (*revIt) {
      case 'A': sequence[i]='T'; break;
      case 'C': sequence[i]='G'; break;
      case 'G': sequence[i]='C'; break;
      case 'T': sequence[i]='A'; break;
      case 'N': sequence[i]='N'; break;
      default: break;
      }
    }
  }
  
  
  template<typename THits, typename TValue>
  inline uint32_t
  findMaxFreq(THits& hits, TValue& gpos) {
    if (hits.empty()) {
      gpos = 0;
      return 0;
    }
    std::sort(hits.begin(), hits.end());
    TValue ihit = hits[0];
    int32_t freq = 1;
    int32_t bestFreq = 1;
    gpos = ihit;
    for(uint32_t i = 1; i<hits.size(); ++i) {
      if (hits[i] == ihit) {
	++freq;
	if (freq > bestFreq) {
	  gpos = ihit;
	  bestFreq = freq;
	}
      } else {
	ihit = hits[i];
	freq = 1;
      }
    }
    return bestFreq;
  }
  

  
  
  template<typename TFMIndex, typename THits>
  inline void
  scanSequence(TFMIndex const& fm_index, std::string const& consensus, uint16_t const trimLeft, uint16_t const trimRight, uint16_t const kmer, THits& hits, bool unique) {
    int32_t ncount = 0;
    for(uint16_t i = trimLeft; ((i < trimLeft + kmer) && (i < consensus.size())); ++i)
      if (consensus[i] == 'N') ++ncount;
    for(uint16_t k = trimLeft + kmer; k < (consensus.size() - trimRight); ++k) {
      if (consensus[k-kmer] == 'N') --ncount;
      if (consensus[k] == 'N') ++ncount;
      if (ncount == 0) {
	std::string seq = consensus.substr(k, kmer);
	std::size_t occs = sdsl::count(fm_index, seq.begin(), seq.end());
	if (unique) {
	  if (occs == 1) {
	    auto locations = locate(fm_index, seq.begin(), seq.end());
	    hits.push_back(locations[0] - k);
	  }
	} else {
	  if (occs > 0) {
	    auto locations = locate(fm_index, seq.begin(), seq.end());
	    for(std::size_t m = 0; m < occs; ++m) {
	      hits.push_back(locations[m] - k);
	    }
	  }
	}
      }
    }
  }


  
  template<typename TConfig, typename TFMIndex>
  inline bool
  getReferenceSlice(TConfig const& c, TFMIndex const& fm_index, teal::BaseCalls const& bc, ReferenceSlice& rs) {
    uint32_t minKmerSupport = 3;
  
    // Get sequence lengths
    std::vector<uint32_t> seqlen;
    faidx_t* fai = NULL;
    if (c.filetype) {
      seqlen.push_back(rs.refslice.size());
    } else {
      fai = fai_load(c.genome.string().c_str());
      seqlen.resize(faidx_nseq(fai));
      for(int32_t refIndex = 0; refIndex < faidx_nseq(fai); ++refIndex) {
	std::string seqname(faidx_iseq(fai, refIndex));
	seqlen[refIndex] = faidx_seq_len(fai, seqname.c_str()) + 1;
      }
    }

    // Fwd and rev index search
    std::vector<int64_t> hitFwd;
    std::vector<int64_t> hitRev;
    scanSequence(fm_index, bc.consensus, c.trimLeft, c.trimRight, c.kmer, hitFwd, true);
    std::string rv = bc.consensus;
    reverseComplement(rv);
    scanSequence(fm_index, rv, c.trimRight, c.trimLeft, c.kmer, hitRev, true);
  
    // Select best orientation
    int64_t bestFwd;
    uint32_t freqFwd = findMaxFreq(hitFwd, bestFwd);
    int64_t bestRev;
    uint32_t freqRev = findMaxFreq(hitRev, bestRev);
    int64_t bestPos;
    if ((freqFwd >= minKmerSupport) && (freqFwd > 2*freqRev)) {
      rs.forward = true;
      rs.kmersupport = freqFwd;
      bestPos = bestFwd;
    } else if ((freqRev >= minKmerSupport) && (freqRev > 2*freqFwd)) {
      rs.forward = false;
      rs.kmersupport = freqRev;
      bestPos = bestRev;
    } else {
      // Try using non-unique matches
      hitFwd.clear();
      hitRev.clear();
      scanSequence(fm_index, bc.consensus, c.trimLeft, c.trimRight, c.kmer, hitFwd, false);
      scanSequence(fm_index, rv, c.trimRight, c.trimLeft, c.kmer, hitRev, false);
      freqFwd = findMaxFreq(hitFwd, bestFwd);
      freqRev = findMaxFreq(hitRev, bestRev);
      if ((freqFwd >= minKmerSupport) && (freqFwd > 2*freqRev)) {
	rs.forward = true;
	rs.kmersupport = freqFwd;
	bestPos = bestFwd;
      } else if ((freqRev >= minKmerSupport) && (freqRev > 2*freqFwd)) {
	rs.forward = false;
	rs.kmersupport = freqRev;
	bestPos = bestRev;
      } else {
	std::cerr << "Couldn't anchor the Sanger trace in the selected reference genome." << std::endl;
	return false;
      }
    }

    // Get initial ref slice
    int64_t cumsum = 0;
    uint32_t refIndex = 0;
    for(; bestPos >= cumsum + seqlen[refIndex]; ++refIndex) cumsum += seqlen[refIndex];
    if (!c.filetype) rs.chr = std::string(faidx_iseq(fai, refIndex));
    uint32_t chrpos = bestPos - cumsum;
    int32_t slen = -1;
    uint32_t slicestart = 0;
    uint32_t sliceend = seqlen[refIndex];
    if (chrpos > c.maxindel) slicestart = chrpos - c.maxindel;
    uint32_t tmpend = chrpos + bc.consensus.size() + c.maxindel;
    if (tmpend < seqlen[refIndex]) sliceend = tmpend;
    if (!c.filetype) {
      rs.pos = slicestart;
      char* seq = faidx_fetch_seq(fai, rs.chr.c_str(), slicestart, sliceend, &slen);
      rs.refslice = boost::to_upper_copy(std::string(seq));
      if (seq != NULL) free(seq);
    }
    if (!rs.forward) reverseComplement(rs.refslice);
    //std::cerr << rs.chr << "\t" << rs.pos << "\t" << rs.forward << std::endl;
    
    // Clean-up
    if (fai != NULL) fai_destroy(fai);
    
    return true;
  }


template<typename TAlign>
inline void
trimReferenceSlice(TAlign const& align, ReferenceSlice& rs) {
  typedef typename TAlign::index TAIndex;
  uint32_t ri = 0;
  int32_t s = -1;
  int32_t e = -1;
  for(TAIndex j = 0; j< (TAIndex) align.shape()[1]; ++j) {
    if (align[0][j] != '-') {
      if (s == -1) s = j;
      e = j + 1;
    }
    if ((s == -1) && (align[1][j] != '-')) ++ri;
  }
  uint32_t risize = 0;
  for(TAIndex j = s; j < (TAIndex) e; ++j) {
    if (align[1][j] != '-') ++risize;
  }
  rs.refslice = rs.refslice.substr(ri, risize);
  rs.pos += ri;
}
  


}

#endif
