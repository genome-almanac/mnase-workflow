/**
 * gaussian.cpp v1.0.0
 *
 * Contains methods related to Gaussian distributions and smoothing.
 * Methods for curve fitting are in gaussianfit.cpp.
 */

#include "gaussian.hpp"

namespace bib {

  namespace gaussian {

    std::vector<double> getDistribution(double A, double x, double u, size_t length) {
      std::vector<double> ret;
      ret.reserve(length + 1);
      for (std::size_t i = 0; i <= length; ++i) {
	ret.push_back(A * std::exp(-(i - x) * (i - x) / (u * u) / 2.0));
      }
      return ret;
    }
  
    /**
     * Samples a Gaussian distribution at a given number of points.
     *
     * stdev: the standard deviation of the distribution.
     * size: the size of the output vector generated will be size * 2 + 1,
     *       with the mean at the center.
     */
    arma::Col<double> sample(double stdev, size_t size) {
      double A = 1.0 / stdev / SQRT2PI;
      arma::Col<double> ret(size * 2 + 1);
      for (auto i = 0; i < size; ++i) {
	ret[size + i] = ret[size - i] = A
	  * std::exp( -i * i / stdev / stdev / 2.0 );
      }
      return ret;
    }

    /**
     * Samples the second derivative of a Gaussian distribution at
     * a given number of points.
     *
     * stdev: the standard deviation of the distribution.
     * size: the size of the output vector generated will be size * 2 + 1,
     *       with the mean at the center.
     */
    arma::Col<double> sampleSDev(double stdev, size_t size) {
      double A = 1.0 / stdev / SQRT2PI;
      arma::Col<double> ret(size * 2 + 1);
      for (auto i = 0; i < size; ++i) {
	ret[size + i] = ret[size - i] = A
	  * std::exp(-i * i / stdev / stdev / 2.0)
	  * (i * i - stdev * stdev) / (stdev * stdev * stdev * stdev);
      }
      return ret;    
    }

    /**
     * Smooths the input vector using a Gaussian distribution of the given width.
     */
    arma::Col<double> smooth(const arma::Col<double>& input, double kernelWidth) {
      auto gaussian = sample(kernelWidth, input.size());
      return arma::conv(input, gaussian, "same");
    }
  
    /**
     * Smooths the input vector using the second derivative of a Gaussian distribution
     * of the given width.
     */
    arma::Col<double> scaleSpaceSmooth(const arma::Col<double>& input, double kernelWidth) {
      auto gaussian = sampleSDev(kernelWidth, input.size());
      return arma::conv(input, gaussian, "same");
    }
    
  }

} // namespace bib::gaussian
