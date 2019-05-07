/**
 * convert.cpp
 * Converts a CSFASTA and CSQUAL pair to FASTQ format.
 */

#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <unordered_map>
#include <vector>

const std::unordered_map<char, std::vector<char>> COLORMAP = {
  { 'A', { 'A', 'C', 'G', 'T' } },
  { 'C', { 'C', 'A', 'T', 'G' } },
  { 'G', { 'G', 'T', 'A', 'C' } },
  { 'T', { 'T', 'G', 'C', 'A' } }
};

std::string decode(const std::string& fastaline) {
  char last = fastaline[0];
  if (COLORMAP.find(last) == COLORMAP.end()) last = 'N';
  std::string out = "";
  for (int i = 1; i < fastaline.size(); ++i) {
    last = fastaline[i] == '.' || last == 'N' ? 'N' : COLORMAP.find(last)->second[fastaline[i] - 48];
    out += last;
  }
  return out;
}

std::string dqual(const std::string& qualline) {
  std::string cint, out = "";
  std::stringstream s(qualline);
  while (std::getline(s, cint, ' ')) {
    out += (char)(cint.compare("-1") == 0 ? 33 : std::atoi(cint.c_str()) + 33);
  }
  return out;
}

int main(int argc, char **argv) {

  if (argc < 4) {
    std::cerr << "usage: fastatofastq input.csfasta input.csqual output.fastq\n";
    return 1;
  }

  size_t i = 0;
  std::ifstream fasta(argv[1]), qual(argv[2]);
  std::string fastaline, qualline, fastaseq, qualseq;
  std::ofstream output(argv[3]);
  std::unordered_map<std::string, std::string> fastas, quals;

  while (std::getline(fasta, fastaline) && std::getline(qual, qualline) && ++i) {
    
    if (i % 100000 == 0) std::cerr << i << '\n';
    
    // skip comments, break on empty line
    while (fastaline[0] == '#' && std::getline(fasta, fastaline));
    while (qualline[0] == '#' && std::getline(qual, qualline));
    if (fastaline.length() == 0 || qualline.length() == 0) break;

    // read actual data
    std::getline(fasta, fastaseq);
    std::getline(qual, qualseq);
    
    // if they are the same, just write to output
    if (fastaline.compare(qualline) == 0) {
      output << '@' << fastaline << '\n' << decode(fastaseq) << "\n+\n" << dqual(qualseq) << '\n';
    } else {

      // if this line from the csqual has already been found in the csfasta, write the pair;
      // otherwise, store it in the csqual map for later
      if (fastas.find(qualline) != fastas.end()) {
	output << '@' << qualline << '\n' << fastas[qualline] << "\n+\n" << dqual(qualseq) << '\n';
	fastas.erase(qualline);
      } else
	quals[qualline] = dqual(qualseq);

      // if this line from the csfasta has already been found in the csqual, write the pair;
      // otherwise, store it in the csfasta map for later
      if (quals.find(fastaline) != quals.end()) {
	output << '@' << fastaline << '\n' << decode(fastaseq) << "\n+\n" << quals[fastaline] << '\n';
	quals.erase(fastaline);
      } else
	fastas[fastaline] = decode(fastaseq);
      
    }
    
  }

  std::cerr << fastaline << '\n' << qualline << '\n';

  return 0;
  
}
