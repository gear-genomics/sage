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

#define BOOST_DISABLE_ASSERTS
#include <boost/multi_array.hpp>
#include <boost/program_options/cmdline.hpp>
#include <boost/program_options/options_description.hpp>
#include <boost/program_options/parsers.hpp>
#include <boost/program_options/variables_map.hpp>
#include <boost/date_time/posix_time/posix_time.hpp>
#include <boost/date_time/gregorian/gregorian.hpp>
#include <boost/iostreams/stream.hpp>
#include <boost/iostreams/stream_buffer.hpp>
#include <boost/iostreams/device/file.hpp>
#include <boost/iostreams/filtering_stream.hpp>
#include <boost/iostreams/filter/zlib.hpp>
#include <boost/iostreams/filter/gzip.hpp>
#include <boost/filesystem.hpp>
#include <boost/progress.hpp>

#include <sdsl/suffix_arrays.hpp>

#include <htslib/faidx.h>

#include "abif.h"
#include "align.h"
#include "gotoh.h"
#include "fmindex.h"
#include "json.h"

using namespace sdsl;
using namespace sage;

struct Config {
  uint16_t trimLeft;
  uint16_t trimRight;
  uint16_t filetype;   //0: *fa.gz, 1: *.fa, 2: *.ab1
  uint16_t kmer;
  uint16_t maxindel;
  uint16_t linelimit;
  float pratio;
  std::string outprefix;
  boost::filesystem::path outfile;
  boost::filesystem::path ab;
  boost::filesystem::path genome;
};

int main(int argc, char** argv) {
  Config c;
  c.outprefix = "bla";
  
  // Parameter
  boost::program_options::options_description generic("Generic options");
  generic.add_options()
    ("help,?", "show help message")
    ("genome,g", boost::program_options::value<boost::filesystem::path>(&c.genome), "(gzipped) fasta or wildtype ab1 file")
    ("pratio,p", boost::program_options::value<float>(&c.pratio)->default_value(0.33), "peak ratio to call base")
    ("kmer,k", boost::program_options::value<uint16_t>(&c.kmer)->default_value(15), "kmer size to anchor trace")
    ("maxindel,m", boost::program_options::value<uint16_t>(&c.maxindel)->default_value(1000), "max. indel size in Sanger trace")
    ("trimLeft,l", boost::program_options::value<uint16_t>(&c.trimLeft)->default_value(50), "trim size left")
    ("trimRight,r", boost::program_options::value<uint16_t>(&c.trimRight)->default_value(50), "trim size right")
    ;

  boost::program_options::options_description otp("Output options");
  otp.add_options()
    ("output,o", boost::program_options::value<boost::filesystem::path>(&c.outfile)->default_value("out.json"), "output file")
    ("linelimit,l", boost::program_options::value<uint16_t>(&c.linelimit)->default_value(60), "alignment line length")
    ;

  boost::program_options::options_description hidden("Hidden options");
  hidden.add_options()
    ("input-file", boost::program_options::value<boost::filesystem::path>(&c.ab), "ab1")
    ;

  boost::program_options::positional_options_description pos_args;
  pos_args.add("input-file", -1);

  boost::program_options::options_description cmdline_options;
  cmdline_options.add(generic).add(otp).add(hidden);
  boost::program_options::options_description visible_options;
  visible_options.add(generic).add(otp);
  boost::program_options::variables_map vm;
  boost::program_options::store(boost::program_options::command_line_parser(argc, argv).options(cmdline_options).positional(pos_args).run(), vm);
  boost::program_options::notify(vm);

  // Check command line arguments
  if ((vm.count("help")) || (!vm.count("input-file"))) {
    std::cout << "Usage: " << argv[0] << " [OPTIONS] -g genome.fa trace.ab1" << std::endl;
    std::cout << visible_options << "\n";
    return -1;
  }
  if (c.maxindel < 1) c.maxindel = 1;

  // Check reference
  if (!(boost::filesystem::exists(c.genome) && boost::filesystem::is_regular_file(c.genome) && boost::filesystem::file_size(c.genome))) {
    std::cerr << "Reference file is missing: " << c.genome.string() << std::endl;
    return 1;
  }

  // Show cmd
  boost::posix_time::ptime now = boost::posix_time::second_clock::local_time();
  std::cout << '[' << boost::posix_time::to_simple_string(now) << "] ";
  for(int i=0; i<argc; ++i) { std::cout << argv[i] << ' '; }
  std::cout << std::endl;

  // Load *.ab1 file
  now = boost::posix_time::second_clock::local_time();
  std::cout << '[' << boost::posix_time::to_simple_string(now) << "] " << "Load ab1 file" << std::endl;
  teal::Trace tr;
  if (!teal::readab(c.ab.string(), tr)) return -1;

  // Call bases
  teal::BaseCalls bc;
  teal::basecall(tr, bc, c.pratio);

  // Load reference
  csa_wt<> fm_index;
  ReferenceSlice rs;
  if (!loadFMIdx(c, rs, fm_index)) return -1;

  // Find reference match
  if (!getReferenceSlice(c, fm_index, bc, rs)) return -1;

  // Semi-global alignment
  typedef boost::multi_array<char, 2> TAlign;
  TAlign align;
  AlignConfig<true, false> semiglobal;
  DnaScore<int> sc(5, -4, -10, -1);
  std::string trace = bc.primary;
  if (trace.size() > (c.trimLeft + c.trimRight)) trace = trace.substr(c.trimLeft, trace.size() - (c.trimLeft + c.trimRight));
  gotoh(trace, rs.refslice, align, semiglobal, sc);

  // Debug Alignment
  //for(uint32_t i = 0; i<align.shape()[0]; ++i) {
  //for(uint32_t j = 0; j<align.shape()[1]; ++j) std::cerr << align[i][j];
  //std::cerr << std::endl;
  //}
  
  // Trim initial reference slice
  trimReferenceSlice(align, rs);

  // Global alignment
  typedef boost::multi_array<char, 2> TAlign;
  TAlign final;
  AlignConfig<false, false> global;
  gotoh(trace, rs.refslice, final, global, sc);

  // Debug Alignment
  //for(uint32_t i = 0; i<final.shape()[0]; ++i) {
  //for(uint32_t j = 0; j<final.shape()[1]; ++j) std::cerr << final[i][j];
  //std::cerr << std::endl;
  //}
  
  // Output
  traceAlignJsonOut(c.outfile.string(), bc, tr, rs, final);
  
  now = boost::posix_time::second_clock::local_time();
  std::cout << '[' << boost::posix_time::to_simple_string(now) << "] " << "Done." << std::endl;

  return 0;
}
