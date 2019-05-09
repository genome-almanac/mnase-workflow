/**
 * bamreader.cpp v1.0.0
 *
 * class for reading and piling up alignments from BAM files.
 */

#include "gaussian.hpp"
#include "bam.hpp"

namespace bib {
  
  namespace bamreader {
    
    BamReader::BamReader(const std::string& path) {
      if (!bamreader_.Open(path)) {
	throw std::invalid_argument("failed to open " + path + " for reading");
      }
      references_ = bamreader_.GetReferenceData();
    }
    
    BamReader::~BamReader() {
      close();
    }
    
    void BamReader::close() {
      bamreader_.Close();
    }

    const std::unordered_map<std::string, std::vector<double>>& BamReader::chromosomeValues() {
      return chromosomes_;
    }
    
    std::unordered_map<std::string, unsigned int> BamReader::chromosomeLengths() {
      std::unordered_map<std::string, unsigned int> ret;
      for (const auto& ref : references_) {
	ret[ref.RefName] = ref.RefLength;
      }
      return ret;
    }
    
    void BamReader::read(char strand) {
    
      BamTools::BamAlignment b;

      for (const auto& ref : references_) {
	chromosomes_[ref.RefName] = std::vector<double>(ref.RefLength);
      }

      const auto ref = references_.at(b.RefID);
      
      if (strand == '.') { /* both strands */
	while (bamreader_.GetNextAlignment(b)) {
	  auto position = b.IsReverseStrand() ? (b.GetEndPosition() - 75) : (b.Position + 75);
	  if (position > 0 && position < chromosomes_[ref.RefName].size()) ++chromosomes_[ref.RefName][position];
	}
      } else if (strand == '-') { /* minus strand only */
	while (bamreader_.GetNextAlignment(b)) {
	  if (b.IsReverseStrand()) {
	    auto position = b.GetEndPosition() - 75;
	    if (position > 0 && position < chromosomes_[ref.RefName].size()) ++chromosomes_[ref.RefName][position];
	  }
	}
      } else if (strand == '+') { /* plus strand only */
	while (bamreader_.GetNextAlignment(b)) {
	  if (!b.IsReverseStrand()) {
	    auto position = b.Position + 75;
	    if (position > 0 && position < chromosomes_[ref.RefName].size()) ++chromosomes_[ref.RefName][position];
	  }
	}
      }
    
    }

    void BamReader::writeWigValues(const arma::Col<double>& values, std::ofstream& f, const std::string& chr) {
      int last_value = -1;
      for (auto i = 0; i < values.size(); ++i) {
	auto& value = values[i];
	if (value <= 0.0) { continue; }
	if (i > last_value) {
	  f << "fixedStep chrom=" << chr << " start=" << (i + 1) << " step=1\n";
	}
	f << value << '\n';
	last_value = i + 1;
      }
    }
    
    void BamReader::writeWig(const std::string& path, double kernel) {
      std::ofstream f(path);
      for (auto& kv : chromosomes_) {
	BamReader::writeWigValues(gaussian::smooth(kv.second, kernel), f, kv.first);
      }
    }
  
  }
  
} // namespace bib::bamreader
