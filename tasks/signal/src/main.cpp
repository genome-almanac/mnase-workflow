#include <unordered_map>
#include <vector>

#include "lib/gaussian.hpp"
#include "lib/bam.hpp"

int main(int argc, char **argv) {

  if (argc < 5) {
    std::cerr << "usage: mnasesignal strand smoothing out.wig *.bam\n";
    return 1;
  }

  double smoothing = std::atof(argv[2]);
  char strand = std::string(argv[1]).compare("+") == 0 ? '+' : (std::string(argv[1]).compare("-") == 0 ? '-' : '.');
  std::unordered_map<std::string, arma::Col<double>> values_;
  for (auto i = 4; i < argc; ++i) {
    bib::bamreader::BamReader b(argv[i]);
    b.read(strand);
    for (auto& kv : b.chromosomeValues()) {
      if (values_.find(kv.first) == values_.end())
	values_[kv.first] = bib::gaussian::smooth(kv.second, smoothing);
      else {
	auto& ovalues = values_[kv.first];
	if (ovalues.size() < kv.second.size())
	  ovalues.resize(kv.second.size());
	auto svalues = bib::gaussian::smooth(kv.second, smoothing);
	for (auto i = 0; i < svalues.size(); ++i)
	  ovalues[i] += svalues[i];
      }
    }
  }

  std::ofstream f(argv[1]);
  for (auto& kv : values_)
    bib::bamreader::BamReader::writeWigValues(kv.second, f, kv.first);
  
  return 0;
  
}
