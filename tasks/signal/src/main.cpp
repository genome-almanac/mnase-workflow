#include <unordered_map>
#include <vector>

#include "lib/gaussian.hpp"
#include "lib/bam.hpp"

int main(int argc, char **argv) {

  if (argc < 3) {
    std::cerr << "usage: mnasesignal out.wig *.bam\n";
    return 1;
  }

  std::unordered_map<std::string, arma::Col<double>> values_;
  for (auto i = 2; i < argc; ++i) {
    bib::bamreader::BamReader b(argv[i]);
    for (auto& kv : b.chromosomeValues()) {
      if (values_.find(kv.first) == values_.end())
	values_[kv.first] = bib::gaussian::smooth(kv.second, 30.0);
      else {
	auto& ovalues = values_[kv.first];
	if (ovalues.size() < kv.second.size())
	  ovalues.resize(kv.second.size());
	auto svalues = bib::gaussian::smooth(kv.second, 30.0);
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
