#ifndef MADGWICK_H_
#define MADGWICK_H_

#include <Eigen/Core>
#include <Eigen/Geometry>

namespace madgwick {

using Vector3 = Eigen::Vector3d;
using Quaternion = Eigen::Quaterniond;

struct DirectionObservation {
  Vector3 reference;
  Vector3 measurement;
};

struct FilterTestAccess;

class Filter {
public:
  explicit Filter(
      double correction_gain = 0.1, double bias_gain = 0.001,
      const Quaternion &initial_orientation = Quaternion::Identity());

  void update(const DirectionObservation &observation,
              const Vector3 &angular_velocity, double dt);

  void update(const DirectionObservation &observation1,
              const DirectionObservation &observation2,
              const Vector3 &angular_velocity, double dt);

  [[nodiscard]] Quaternion orientation() const noexcept;
  void set_orientation(const Quaternion &orientation);

  [[nodiscard]] Vector3 gyro_bias() const noexcept;
  void set_gyro_bias(const Vector3 &bias);

  [[nodiscard]] double correction_gain() const noexcept;
  [[nodiscard]] double bias_gain() const noexcept;
  void set_gains(double correction_gain, double bias_gain);

  void reset(const Quaternion &orientation = Quaternion::Identity(),
             const Vector3 &gyro_bias = Vector3::Zero());

  [[nodiscard]] static Vector3 angular_velocity(const Quaternion &current,
                                                const Quaternion &previous,
                                                double dt);

private:
  friend struct FilterTestAccess;

  [[nodiscard]] static Quaternion
  observation_gradient(const DirectionObservation &observation,
                       const Quaternion &orientation);

  void apply_gradient(const Quaternion &gradient,
                      const Vector3 &angular_velocity, double dt);

  double correction_gain_;
  double bias_gain_;
  Quaternion orientation_ = Quaternion::Identity();
  Vector3 gyro_bias_ = Vector3::Zero();
};

} // namespace madgwick

#endif
