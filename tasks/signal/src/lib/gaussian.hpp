/**
 * gaussian.hpp v1.0.0
 *
 * Contains methods related to Gaussian distributions and smoothing.
 * Methods for curve fitting are in gaussianfit.cpp.
 */

#pragma once

#define SQRT2PI 2.506628275

#include <cmath>

#define  ARMA_DONT_USE_WRAPPER
#define  ARMA_USE_LAPACK
#define  ARMA_DONT_PRINT_ERRORS
#include <armadillo>

namespace bib {

  namespace gaussian {

    const double M_SQRT_2OVERPI = std::sqrt(2.0 / M_PI);
    const double M_SQRT_2OVERPI_CUBED = std::pow(M_SQRT_2OVERPI, 3.0);
  
    /**
     * Represents a Gaussian distribution. Comparison functions only
     * compare amplitude (used for sorting from most to least important
     * when curve fitting).
     *
     * amplitude: suqare root of the amplitude of the distribution.
     * mean: center point of the distribution.
     * stdev: standard deviation of the distribution.
     */
    struct Parameters {
    
      double amplitude;
      double mean;
      double stdev;

      bool operator >(const Parameters& a) const {
	return amplitude * amplitude / stdev > a.amplitude * a.amplitude / a.stdev;
      }
    
      bool operator <(const Parameters& a) const {
	return amplitude * amplitude / stdev < a.amplitude * a.amplitude / a.stdev;
      }
    
    };

    /**
     * Represents a Gaussian distribution with shape (adds skewness).
     * Comparison functions only compare amplitude (used for sorting
     * from most to least important when curve fitting).
     *
     * shape: shape parameter (alpha) capturing skewness.
     * amplitude: suqare root of the amplitude of the distribution.
     * mean: center point of the distribution.
     * stdev: standard deviation of the distribution.
     */
    struct SkewParameters {
    
      double amplitude;
      double mean;
      double stdev;
      double shape;

      bool operator >(const SkewParameters& a) const {
	return amplitude * amplitude / stdev > a.amplitude * a.amplitude / a.stdev;
      }
    
      bool operator <(const SkewParameters& a) const {
	return amplitude * amplitude / stdev < a.amplitude * a.amplitude / a.stdev;
      }
    
    };

    /**
     * Returns the value of a Gaussian distribution at a given point x.
     *
     * amplitude: suqare root of the amplitude of the distribution.
     * mean: center point of the distribution.
     * stdev: standard deviation of the distribution.
     * x: the point at which to sample the distribution.
     */
    inline double curveValue(double amplitude, double mean, double stdev, double x) {
      return amplitude * amplitude / stdev * std::exp( -(x - mean) * (x - mean) / stdev / stdev / 2.0 );
    }

    /**
     * Returns the value of a Gaussian distribution with shape at a given point x.
     *
     * amplitude: suqare root of the amplitude of the distribution.
     * mean: center point of the distribution.
     * stdev: standard deviation of the distribution.
     * shape: shape parameter (alpha) capturing skewness.
     * x: the point at which to sample the distribution.
     */
    inline double curveValue(double amplitude, double mean, double stdev, double shape, double x) {
      return amplitude * amplitude / stdev * std::exp( -(x - mean) * (x - mean) / stdev / stdev / 2.0 )
	* (1.0 + std::erf(shape * (x - mean) / (stdev * M_SQRT2)));
    }
  
    /**
     * Returns the value of a Gaussian distribution with shape at a given point x.
     *
     * distribution: the distribution to sample.
     * x: the point at which to sample the distribution.
     */
    inline double curveValue(Parameters distribution, double x) {
      return curveValue(distribution.amplitude, distribution.mean, distribution.stdev, x);
    }

    /**
     * Returns the value of a Gaussian distribution with shape at a given point x.
     *
     * distribution: the distribution to sample.
     * x: the point at which to sample the distribution.
     */
    inline double curveValue(SkewParameters distribution, double x) {
      return curveValue(distribution.amplitude, distribution.mean, distribution.stdev,
			distribution.shape, x);
    }
  
    /**
     * Smooths the input vector using a Gaussian distribution of the given width.
     *
     * input: the input vector to smooth.
     * kernelWidth: standard deviation of the distribution to use for smoothing.
     */
    arma::Col<double> smooth(const arma::Col<double>& input, double kernelWidth);

    /**
     * Smooths the input vector using the second derivative of a Gaussian distribution
     * of the given width.
     *
     * input: the input vector to smooth.
     * kernelWidth: standard deviation of the distribution to use for smoothing.
     */
    arma::Col<double> scaleSpaceSmooth(const arma::Col<double>& input, double kernelWidth);

    std::vector<double> getDistribution(double A, double x, double u, size_t length);

    inline int sgn(double val) {
      return (0.0 < val) - (val < 0.0);
    }
  
    inline double skewMode(SkewParameters gaussian) {
      double delta = gaussian.shape / std::sqrt(1.0 + gaussian.shape * gaussian.shape);
      double uz = M_SQRT_2OVERPI * delta;
      double oz = std::sqrt(1.0 - uz * uz);
      double skewness = ((4.0 - M_PI) / 2.0)
	* ((M_SQRT_2OVERPI_CUBED * delta * delta * delta) / std::pow(1.0 - 2.0 * delta * delta / M_PI, 1.5));
      return (uz - skewness * oz / 2.0 - sgn(gaussian.shape) / 2.0 * std::exp(-2.0 * M_PI / std::abs(gaussian.shape)))
	* gaussian.stdev + gaussian.mean;
    }
    
  }

} // namespace bib::gaussian
