#include "Madgwick.h"

#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <string>

namespace madgwick {
namespace {

constexpr double minimum_norm = 64.0 * std::numeric_limits<double>::epsilon();

Vector3 normalized_direction(const Vector3 &value, const char *name) {
  const double norm = value.stableNorm();
  if (!value.allFinite() || !std::isfinite(norm) || norm <= minimum_norm)
    throw std::invalid_argument(std::string(name) +
                                " must be finite and non-zero");
  return value / norm;
}

Quaternion normalized_quaternion(const Quaternion &value, const char *name) {
  const double norm = value.coeffs().stableNorm();
  if (!value.coeffs().allFinite() || !std::isfinite(norm) ||
      norm <= minimum_norm)
    throw std::invalid_argument(std::string(name) +
                                " must be finite and non-zero");
  Quaternion result = value;
  result.coeffs() /= norm;
  return result;
}

void validate_gain(double value, const char *name) {
  if (!std::isfinite(value) || value < 0.0)
    throw std::invalid_argument(std::string(name) +
                                " must be finite and non-negative");
}

void validate_step(const Vector3 &angular_velocity, double dt) {
  if (!angular_velocity.allFinite())
    throw std::invalid_argument("angular velocity must be finite");
  if (!std::isfinite(dt) || dt <= 0.0)
    throw std::invalid_argument("dt must be finite and positive");
}

} // namespace

Filter::Filter(double correction_gain, double bias_gain,
               const Quaternion &initial_orientation)
    : correction_gain_(correction_gain), bias_gain_(bias_gain) {
  set_gains(correction_gain, bias_gain);
  set_orientation(initial_orientation);
}

void Filter::update(const DirectionObservation &observation,
                    const Vector3 &angular_velocity, double dt) {
  validate_step(angular_velocity, dt);
  const DirectionObservation normalized_observation{
      normalized_direction(observation.reference, "reference"),
      normalized_direction(observation.measurement, "measurement")};
  apply_gradient(observation_gradient(normalized_observation, orientation_),
                 angular_velocity, dt);
}

void Filter::update(const DirectionObservation &observation1,
                    const DirectionObservation &observation2,
                    const Vector3 &angular_velocity, double dt) {
  validate_step(angular_velocity, dt);
  const DirectionObservation normalized_observation1{
      normalized_direction(observation1.reference, "reference1"),
      normalized_direction(observation1.measurement, "measurement1")};
  const DirectionObservation normalized_observation2{
      normalized_direction(observation2.reference, "reference2"),
      normalized_direction(observation2.measurement, "measurement2")};

  Quaternion gradient =
      observation_gradient(normalized_observation1, orientation_);
  gradient.coeffs() +=
      observation_gradient(normalized_observation2, orientation_).coeffs();
  apply_gradient(gradient, angular_velocity, dt);
}

Quaternion Filter::orientation() const noexcept {
  return orientation_.conjugate();
}

void Filter::set_orientation(const Quaternion &orientation) {
  orientation_ = normalized_quaternion(orientation, "orientation").conjugate();
}

Vector3 Filter::gyro_bias() const noexcept { return gyro_bias_; }

void Filter::set_gyro_bias(const Vector3 &bias) {
  if (!bias.allFinite())
    throw std::invalid_argument("gyroscope bias must be finite");
  gyro_bias_ = bias;
}

double Filter::correction_gain() const noexcept { return correction_gain_; }

double Filter::bias_gain() const noexcept { return bias_gain_; }

void Filter::set_gains(double correction_gain, double bias_gain) {
  validate_gain(correction_gain, "correction gain");
  validate_gain(bias_gain, "bias gain");
  correction_gain_ = correction_gain;
  bias_gain_ = bias_gain;
}

void Filter::reset(const Quaternion &orientation, const Vector3 &gyro_bias) {
  const Quaternion normalized_orientation =
      normalized_quaternion(orientation, "orientation");
  if (!gyro_bias.allFinite())
    throw std::invalid_argument("gyroscope bias must be finite");
  orientation_ = normalized_orientation.conjugate();
  gyro_bias_ = gyro_bias;
}

Vector3 Filter::angular_velocity(const Quaternion &current,
                                 const Quaternion &previous, double dt) {
  if (!std::isfinite(dt) || dt <= 0.0)
    throw std::invalid_argument("dt must be finite and positive");

  const Quaternion normalized_current =
      normalized_quaternion(current, "current orientation");
  const Quaternion normalized_previous =
      normalized_quaternion(previous, "previous orientation");
  Quaternion delta = normalized_previous.conjugate() * normalized_current;
  delta.normalize();

  if (delta.w() < 0.0)
    delta.coeffs() *= -1.0;

  const double vector_norm = delta.vec().norm();
  Vector3 result;
  if (vector_norm <= minimum_norm) {
    result = 2.0 * delta.vec() / dt;
  } else {
    const double angle =
        2.0 * std::atan2(vector_norm, std::clamp(delta.w(), 0.0, 1.0));
    result = angle * delta.vec() / (vector_norm * dt);
  }

  if (!result.allFinite())
    throw std::overflow_error("angular velocity is not representable");
  return result;
}

Quaternion Filter::observation_gradient(const DirectionObservation &observation,
                                        const Quaternion &orientation) {
  const double w = orientation.w();
  const double x = orientation.x();
  const double y = orientation.y();
  const double z = orientation.z();
  const double dx = observation.reference.x();
  const double dy = observation.reference.y();
  const double dz = observation.reference.z();

  const Vector3 predicted = orientation * observation.reference;
  const Vector3 residual = predicted - observation.measurement;

  Eigen::Matrix<double, 3, 4> jacobian;
  jacobian << 2.0 * (w * dx - z * dy + y * dz),
      2.0 * (x * dx + y * dy + z * dz), 2.0 * (-y * dx + x * dy + w * dz),
      2.0 * (-z * dx - w * dy + x * dz), 2.0 * (z * dx + w * dy - x * dz),
      2.0 * (y * dx - x * dy - w * dz), 2.0 * (x * dx + y * dy + z * dz),
      2.0 * (w * dx - z * dy + y * dz), 2.0 * (-y * dx + x * dy + w * dz),
      2.0 * (z * dx + w * dy - x * dz), 2.0 * (-w * dx + z * dy - y * dz),
      2.0 * (x * dx + y * dy + z * dz);

  const Eigen::Vector4d coefficients = jacobian.transpose() * residual;
  return Quaternion(coefficients[0], coefficients[1], coefficients[2],
                    coefficients[3]);
}

void Filter::apply_gradient(const Quaternion &gradient,
                            const Vector3 &angular_velocity, double dt) {
  const Quaternion bias_error_quaternion = gradient * orientation_.conjugate();
  const Vector3 next_bias =
      gyro_bias_ - 2.0 * bias_gain_ * bias_error_quaternion.vec() * dt;

  Quaternion angular_velocity_quaternion(0.0,
                                         angular_velocity.x() - next_bias.x(),
                                         angular_velocity.y() - next_bias.y(),
                                         angular_velocity.z() - next_bias.z());
  Quaternion derivative = angular_velocity_quaternion * orientation_;
  derivative.coeffs() *= -0.5;
  derivative.coeffs() -= correction_gain_ * gradient.coeffs();

  Quaternion next_orientation = orientation_;
  next_orientation.coeffs() += derivative.coeffs() * dt;
  const double next_orientation_norm = next_orientation.coeffs().stableNorm();
  if (!next_bias.allFinite() || !next_orientation.coeffs().allFinite() ||
      !std::isfinite(next_orientation_norm) ||
      next_orientation_norm <= minimum_norm)
    throw std::overflow_error("filter update produced an invalid state");

  next_orientation.coeffs() /= next_orientation_norm;
  gyro_bias_ = next_bias;
  orientation_ = next_orientation;
}

} // namespace madgwick
