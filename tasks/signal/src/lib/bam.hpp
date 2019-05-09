/**
 * bamreader.hpp v1.0.0
 *
 * class for reading and piling up alignments from BAM files.
 */

#pragma once

#include <unordered_map>
#include <stdexcept>

#include <bamtools/api/BamAlignment.h>
#include <bamtools/api/BamReader.h>

namespace bib {

  namespace bamreader {

    /**
     * Provides methods for reading BAM files and piling up alignments or
     * their 5' cut sites.
     *
     * bamreader_: the underlying bam reader object from BamTools.
     * references_: vector of reference sequence (chromosome) objects present in the bam.
     * values_: after read or read_atac is called, contains alignment pileup positions.
     */
    class BamReader {

    private:
      std::unordered_map<std::string, std::vector<double>> chromosomes_;
    
    protected:
      BamTools::BamReader bamreader_;
      BamTools::RefVector references_;
    
    public:

      /**
       * Constructs a new BamReader object; opens the BAM at the given path and
       * reads the list of reference sequences. The BAM remains open for reading
       * alignments and is closed by close() or the descructor.
       *
       * path: path to the input BAM.
       */
      BamReader(const std::string& path);

      void writeWig(const std::string& path, double kernel);

      static void writeWigValues(const arma::Col<double>& values, std::ofstream& f, const std::string& chr);

      const std::unordered_map<std::string, std::vector<double>>& chromosomeValues();

      /**
       * Destructor; closes the open BAM file.
       */
      virtual ~BamReader();

      /**
       * Closes the open BAM file; called by the destructor.
       */
      void close();

      /**
       * Returns a map of chromosome names to their lengths.
       */
      std::unordered_map<std::string, unsigned int> chromosomeLengths();

      /**
       * Returns a map of chromosome names to the coordinates of the piled up
       * alignments of 5' cut locations.
       */
      const std::unordered_map<std::string, std::vector<unsigned int>>& values();

      /**
       * Reads all alignments in the BAM file and piles up the 5' cut sites.
       *
       * strand: if '+' or '-', only reads alignments from the plus or minus strand.
       */
      void read(char strand = '.');
    
    };

  }

} // namespace bib::bamreader
