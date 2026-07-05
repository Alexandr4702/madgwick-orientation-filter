#include "Madgwick.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <gtest/gtest.h>
#include <limits>
#include <random>
#include <stdexcept>

using madgwick::DirectionObservation;
using madgwick::Filter;
using madgwick::Quaternion;
using madgwick::Vector3;

namespace madgwick {
struct FilterTestAccess {
  static Quaternion observationGradient(const Vector3 &reference,
                                        const Vector3 &measurement,
                                        const Quaternion &orientation) {
    return Filter::observation_gradient({reference, measurement}, orientation);
  }
};
} // namespace madgwick

namespace {
DirectionObservation observation(const Vector3 &reference,
                                 const Vector3 &measurement) {
  return {reference, measurement};
}

double angularError(const Quaternion &actual, const Quaternion &expected) {
  const double dot =
      std::min(1.0, std::abs(actual.normalized().dot(expected.normalized())));
  return 2.0 * std::acos(dot);
}

Quaternion integrateTruth(const Quaternion &orientation, const Vector3 &omega,
                          double dt) {
  const double angle = omega.norm() * dt;
  if (angle == 0.0)
    return orientation;
  return (orientation *
          Quaternion(Eigen::AngleAxisd(angle, omega.normalized())))
      .normalized();
}

Vector3 noisyDirection(const Vector3 &direction, std::mt19937 &generator,
                       std::normal_distribution<double> &noise) {
  return (direction +
          Vector3(noise(generator), noise(generator), noise(generator)))
      .normalized();
}

double observationCost(const Quaternion &orientation, const Vector3 &reference,
                       const Vector3 &measurement) {
  const Quaternion direction(0.0, reference.x(), reference.y(), reference.z());
  const Vector3 residual =
      (orientation * direction * orientation.conjugate()).vec() - measurement;
  return 0.5 * residual.squaredNorm();
}

Quaternion perturbed(const Quaternion &orientation, int component,
                     double amount) {
  Quaternion result = orientation;
  if (component == 0)
    result.w() += amount;
  else
    result.coeffs()[component - 1] += amount;
  return result;
}
} // namespace

TEST(Characterization, IdentityIsPreservedWithZeroAngularVelocity) {
  Filter filter;
  for (int i = 0; i < 100; ++i)
    filter.update(observation(Vector3::UnitX(), Vector3::UnitX()),
                  observation(Vector3::UnitY(), Vector3::UnitY()),
                  Vector3::Zero(), 0.01);

  const Quaternion q = filter.orientation();
  EXPECT_DOUBLE_EQ(q.w(), 1.0);
  EXPECT_DOUBLE_EQ(q.x(), 0.0);
  EXPECT_DOUBLE_EQ(q.y(), 0.0);
  EXPECT_DOUBLE_EQ(q.z(), 0.0);
}

TEST(Characterization, OneGyroStepMatchesExistingEulerIntegration) {
  Filter filter;
  filter.update(observation(Vector3::UnitZ(), Vector3::UnitZ()),
                Vector3(1.0, 2.0, 3.0), 0.1);

  // q_BW integrates (1, -0.05, -0.10, -0.15); orientation() returns q_BW^-1.
  const double norm = std::sqrt(1.0 + 0.05 * 0.05 + 0.10 * 0.10 + 0.15 * 0.15);
  const Quaternion q = filter.orientation();
  EXPECT_NEAR(q.w(), 1.0 / norm, 1e-15);
  EXPECT_NEAR(q.x(), 0.05 / norm, 1e-15);
  EXPECT_NEAR(q.y(), 0.10 / norm, 1e-15);
  EXPECT_NEAR(q.z(), 0.15 / norm, 1e-15);
}

TEST(Kinematics, NonIdentityPropagationMatchesBodyFrameAngularVelocity) {
  const Quaternion initial(
      Eigen::AngleAxisd(0.8, Vector3(1.0, 2.0, -1.0).normalized()));
  const Vector3 omega(0.3, -0.2, 0.1);
  constexpr double dt = 1e-4;
  Filter filter(0.0, 0.0, initial);

  filter.update(
      observation(Vector3::UnitZ(), initial.conjugate() * Vector3::UnitZ()),
      omega, dt);

  const Quaternion expected =
      (initial *
       Quaternion(Eigen::AngleAxisd(omega.norm() * dt, omega.normalized())))
          .normalized();
  const Quaternion actual = filter.orientation();
  const Vector3 recovered = Filter::angular_velocity(actual, initial, dt);

  EXPECT_LT(angularError(actual, expected), 1e-7);
  EXPECT_LT((recovered - omega).norm(), 1e-8);
}

TEST(Math, AnalyticObservationGradientMatchesFiniteDifferences) {
  std::mt19937 generator(0x47524144);
  std::normal_distribution<double> distribution;
  constexpr double epsilon = 1e-7;

  for (int sample = 0; sample < 100; ++sample) {
    const Quaternion orientation(
        distribution(generator), distribution(generator),
        distribution(generator), distribution(generator));
    const Quaternion normalized_orientation = orientation.normalized();
    const Vector3 reference =
        Vector3(distribution(generator), distribution(generator),
                distribution(generator))
            .normalized();
    const Vector3 measurement =
        Vector3(distribution(generator), distribution(generator),
                distribution(generator))
            .normalized();
    const Quaternion analytic = madgwick::FilterTestAccess::observationGradient(
        reference, measurement, normalized_orientation);
    const std::array<double, 4> analytic_components = {
        analytic.w(), analytic.x(), analytic.y(), analytic.z()};

    for (int component = 0; component < 4; ++component) {
      const double numeric = (observationCost(perturbed(normalized_orientation,
                                                        component, epsilon),
                                              reference, measurement) -
                              observationCost(perturbed(normalized_orientation,
                                                        component, -epsilon),
                                              reference, measurement)) /
                             (2.0 * epsilon);
      EXPECT_NEAR(analytic_components[component], numeric, 2e-8)
          << "sample " << sample << ", component " << component;
    }
  }
}

TEST(Simulation, TracksConstantRotationWithIdealSensors) {
  Filter filter;
  Quaternion truth = Quaternion::Identity();
  const Vector3 omega(0.08, -0.12, 0.20);
  const Vector3 ref1 = Vector3::UnitZ();
  const Vector3 ref2 = Vector3(0.8, 0.2, 0.1).normalized();
  constexpr double dt = 0.002;

  for (int i = 0; i < 5000; ++i) {
    truth = integrateTruth(truth, omega, dt);
    const Vector3 meas1 = truth.conjugate() * ref1;
    const Vector3 meas2 = truth.conjugate() * ref2;
    filter.update(observation(ref1, meas1), observation(ref2, meas2), omega,
                  dt);
  }

  EXPECT_LT(angularError(filter.orientation(), truth), 5e-4);
}

TEST(Simulation, OneVectorUpdateTracksObservableMotion) {
  Filter filter;
  Quaternion truth = Quaternion::Identity();
  const Vector3 reference = Vector3::UnitZ();
  const Vector3 omega(0.2, -0.1, 0.0);
  constexpr double dt = 0.001;

  for (int i = 0; i < 10000; ++i) {
    truth = integrateTruth(truth, omega, dt);
    filter.update(observation(reference, truth.conjugate() * reference), omega,
                  dt);
  }

  EXPECT_LT(angularError(filter.orientation(), truth), 1e-3);
}

TEST(Simulation, CorrectsAnInitialAttitudeError) {
  Filter filter;
  const Quaternion truth(
      Eigen::AngleAxisd(0.35, Vector3(1.0, -2.0, 0.5).normalized()));
  const Vector3 ref1 = Vector3::UnitZ();
  const Vector3 ref2 = Vector3::UnitX();
  const Vector3 meas1 = truth.conjugate() * ref1;
  const Vector3 meas2 = truth.conjugate() * ref2;

  for (int i = 0; i < 5000; ++i)
    filter.update(observation(ref1, meas1), observation(ref2, meas2),
                  Vector3::Zero(), 0.002);

  // The unnormalized gradient settles to a small residual error.
  EXPECT_LT(angularError(filter.orientation(), truth), 5e-3);
}

TEST(Simulation, EstimatesConstantGyroscopeBiasWhileStationary) {
  Filter filter(0.1, 0.01);
  const Vector3 gyro_bias(0.004, -0.003, 0.002);
  const Vector3 ref1 = Vector3::UnitZ();
  const Vector3 ref2 = Vector3::UnitX();

  for (int i = 0; i < 30000; ++i)
    filter.update(observation(ref1, ref1), observation(ref2, ref2), gyro_bias,
                  0.002);

  EXPECT_LT(angularError(filter.orientation(), Quaternion::Identity()), 1e-3);
  EXPECT_NEAR(filter.gyro_bias().x(), gyro_bias.x(), 5e-4);
  EXPECT_NEAR(filter.gyro_bias().y(), gyro_bias.y(), 5e-4);
  EXPECT_NEAR(filter.gyro_bias().z(), gyro_bias.z(), 5e-4);
}

TEST(Simulation, GyroscopeBiasConvergesFromAnIncorrectInitialEstimate) {
  Filter filter(0.1, 0.01);
  const Vector3 actual_bias(-0.005, 0.003, 0.004);
  filter.set_gyro_bias(Vector3(0.003, -0.002, 0.001));
  const Vector3 ref1 = Vector3(0.2, -0.3, 0.93).normalized();
  const Vector3 ref2 = Vector3(-0.7, 0.6, 0.2).normalized();

  for (int i = 0; i < 40000; ++i)
    filter.update(observation(ref1, ref1), observation(ref2, ref2), actual_bias,
                  0.002);

  EXPECT_LT((filter.gyro_bias() - actual_bias).norm(), 5e-4);
  EXPECT_LT(angularError(filter.orientation(), Quaternion::Identity()), 1e-3);
}

TEST(Simulation, ZeroBiasGainLeavesTheBiasEstimateUnchanged) {
  const Vector3 initial_bias(0.003, -0.002, 0.001);
  Filter filter(0.1, 0.0);
  filter.set_gyro_bias(initial_bias);
  const Vector3 ref1 = Vector3::UnitZ();
  const Vector3 ref2 = Vector3::UnitX();

  for (int i = 0; i < 1000; ++i)
    filter.update(observation(ref1, ref1), observation(ref2, ref2),
                  Vector3(0.01, -0.02, 0.005), 0.002);

  EXPECT_EQ(filter.gyro_bias(), initial_bias);
}

TEST(Simulation, TracksTimeVaryingThreeAxisMotionWithSamplePeriodJitter) {
  Filter filter;
  Quaternion truth = Quaternion::Identity();
  const Vector3 ref1 = Vector3(0.2, -0.3, 0.93).normalized();
  const Vector3 ref2 = Vector3(-0.7, 0.6, 0.2).normalized();
  double time = 0.0;
  double maximum_error = 0.0;

  for (int i = 0; i < 12000; ++i) {
    const double dt = 0.0015 + 0.0005 * (1.0 + std::sin(0.017 * i));
    const double middle = time + 0.5 * dt;
    const Vector3 omega(0.45 * std::sin(0.71 * middle),
                        -0.35 * std::cos(0.43 * middle),
                        0.25 + 0.20 * std::sin(0.19 * middle));

    truth = integrateTruth(truth, omega, dt);
    filter.update(observation(ref1, truth.conjugate() * ref1),
                  observation(ref2, truth.conjugate() * ref2), omega, dt);
    maximum_error =
        std::max(maximum_error, angularError(filter.orientation(), truth));
    time += dt;
  }

  EXPECT_LT(angularError(filter.orientation(), truth), 2e-3);
  EXPECT_LT(maximum_error, 3e-3);
}

TEST(Simulation, RemainsAccurateWithGyroscopeAndDirectionNoise) {
  Filter filter;
  Quaternion truth = Quaternion::Identity();
  const Vector3 ref1 = Vector3::UnitZ();
  const Vector3 ref2 = Vector3(0.8, -0.3, 0.15).normalized();
  std::mt19937 generator(0x4d414447);
  std::normal_distribution<double> direction_noise(0.0, 0.004);
  std::normal_distribution<double> gyro_noise(0.0, 0.002);
  constexpr double dt = 0.002;
  double maximum_error = 0.0;

  for (int i = 0; i < 15000; ++i) {
    const double time = i * dt;
    const Vector3 omega(0.18 * std::sin(0.37 * time),
                        0.12 * std::cos(0.23 * time),
                        -0.16 * std::sin(0.51 * time));
    truth = integrateTruth(truth, omega, dt);

    const Vector3 gyro =
        omega + Vector3(gyro_noise(generator), gyro_noise(generator),
                        gyro_noise(generator));
    const Vector3 meas1 =
        noisyDirection(truth.conjugate() * ref1, generator, direction_noise);
    const Vector3 meas2 =
        noisyDirection(truth.conjugate() * ref2, generator, direction_noise);
    filter.update(observation(ref1, meas1), observation(ref2, meas2), gyro, dt);
    maximum_error =
        std::max(maximum_error, angularError(filter.orientation(), truth));
  }

  EXPECT_LT(angularError(filter.orientation(), truth), 1e-2);
  EXPECT_LT(maximum_error, 2e-2);
}

TEST(Simulation, TracksManeuversWhileEstimatingConstantGyroscopeBias) {
  Filter filter;
  Quaternion truth = Quaternion::Identity();
  const Vector3 ref1 = Vector3(0.3, 0.1, 0.95).normalized();
  const Vector3 ref2 = Vector3(-0.2, 0.9, 0.35).normalized();
  const Vector3 gyro_bias(0.006, -0.004, 0.003);
  constexpr double dt = 0.002;
  double maximum_error = 0.0;

  for (int i = 0; i < 30000; ++i) {
    const double time = i * dt;
    const Vector3 omega(0.20 * std::sin(0.29 * time),
                        -0.15 * std::sin(0.47 * time),
                        0.10 * std::cos(0.31 * time));
    truth = integrateTruth(truth, omega, dt);
    filter.update(observation(ref1, truth.conjugate() * ref1),
                  observation(ref2, truth.conjugate() * ref2),
                  omega + gyro_bias, dt);
    maximum_error =
        std::max(maximum_error, angularError(filter.orientation(), truth));
  }

  EXPECT_LT(angularError(filter.orientation(), truth), 2e-2);
  EXPECT_LT(maximum_error, 3e-2);
  EXPECT_LT((filter.gyro_bias() - gyro_bias).norm(), gyro_bias.norm());
}

TEST(Simulation, ConvergesFromLargeErrorsAroundDifferentAxes) {
  const Vector3 ref1 = Vector3(0.1, -0.2, 0.97).normalized();
  const Vector3 ref2 = Vector3(0.8, 0.5, -0.1).normalized();
  const std::array<Vector3, 4> axes = {Vector3::UnitX(), Vector3::UnitY(),
                                       Vector3::UnitZ(),
                                       Vector3(1.0, -2.0, 3.0).normalized()};

  for (const Vector3 &axis : axes) {
    Filter filter;
    const Quaternion truth(Eigen::AngleAxisd(1.4, axis));
    const Vector3 meas1 = truth.conjugate() * ref1;
    const Vector3 meas2 = truth.conjugate() * ref2;

    for (int i = 0; i < 20000; ++i)
      filter.update(observation(ref1, meas1), observation(ref2, meas2),
                    Vector3::Zero(), 0.002);

    EXPECT_LT(angularError(filter.orientation(), truth), 3e-2)
        << "axis: " << axis.transpose();
  }
}

TEST(Characterization, QuaternionDeltaProducesExpectedAngularVelocity) {
  const double angle = 0.2;
  const Quaternion previous = Quaternion::Identity();
  const Quaternion current(Eigen::AngleAxisd(angle, Vector3::UnitY()));
  const Vector3 omega = Filter::angular_velocity(current, previous, 0.1);
  EXPECT_NEAR(omega.x(), 0.0, 1e-14);
  EXPECT_NEAR(omega.y(), 2.0, 1e-13);
  EXPECT_NEAR(omega.z(), 0.0, 1e-14);
}

TEST(Characterization, QuaternionDeltaUsesShortestEquivalentRotation) {
  Quaternion current(Eigen::AngleAxisd(1e-6, Vector3::UnitX()));
  current.coeffs() *= -1.0;

  const Vector3 omega =
      Filter::angular_velocity(current, Quaternion::Identity(), 0.01);

  EXPECT_NEAR(omega.x(), 1e-4, 1e-12);
  EXPECT_NEAR(omega.y(), 0.0, 1e-14);
  EXPECT_NEAR(omega.z(), 0.0, 1e-14);
}

TEST(Validation, RejectsInvalidInputsWithoutChangingState) {
  Filter filter;
  const Quaternion initial_orientation = filter.orientation();
  const Vector3 initial_bias = filter.gyro_bias();
  const double nan = std::numeric_limits<double>::quiet_NaN();

  EXPECT_THROW(filter.update(observation(Vector3::Zero(), Vector3::UnitZ()),
                             Vector3::Zero(), 0.01),
               std::invalid_argument);
  EXPECT_THROW(filter.update(observation(Vector3::UnitZ(), Vector3::UnitZ()),
                             Vector3::Zero(), 0.0),
               std::invalid_argument);
  EXPECT_THROW(filter.update(observation(Vector3::UnitZ(), Vector3::UnitZ()),
                             Vector3(nan, 0.0, 0.0), 0.01),
               std::invalid_argument);
  EXPECT_THROW((void)Filter::angular_velocity(Quaternion::Identity(),
                                              Quaternion::Identity(), -0.01),
               std::invalid_argument);

  EXPECT_EQ(filter.orientation().coeffs(), initial_orientation.coeffs());
  EXPECT_EQ(filter.gyro_bias(), initial_bias);
}

TEST(Validation, NormalizesDirectionsAndInitialOrientation) {
  const Quaternion initial(2.0, 0.0, 0.0, 0.0);
  Filter filter(0.1, 0.001, initial);
  filter.update(observation(5.0 * Vector3::UnitZ(), 3.0 * Vector3::UnitZ()),
                Vector3::Zero(), 0.01);

  EXPECT_NEAR(filter.orientation().norm(), 1.0, 1e-15);
  EXPECT_LT(angularError(filter.orientation(), Quaternion::Identity()), 1e-15);
}

TEST(Validation, HandlesLargeFiniteDirectionMagnitudes) {
  Filter filter;
  const double large = std::numeric_limits<double>::max() / 4.0;

  EXPECT_NO_THROW(filter.update(
      observation(Vector3(large, 0.0, 0.0), Vector3(large, 0.0, 0.0)),
      Vector3::Zero(), 0.01));
  EXPECT_TRUE(filter.orientation().coeffs().allFinite());
}

TEST(State, BiasGainsAndResetRoundTrip) {
  Filter filter(0.2, 0.003);
  const Vector3 bias(0.125, -0.25, 0.5);
  filter.set_gyro_bias(bias);
  const Quaternion orientation(Eigen::AngleAxisd(0.4, Vector3::UnitY()));
  filter.set_orientation(orientation);

  EXPECT_EQ(filter.gyro_bias(), bias);
  EXPECT_NEAR(filter.correction_gain(), 0.2, 0.0);
  EXPECT_NEAR(filter.bias_gain(), 0.003, 0.0);
  EXPECT_LT(angularError(filter.orientation(), orientation), 1e-15);

  filter.reset();
  EXPECT_EQ(filter.gyro_bias(), Vector3::Zero());
  EXPECT_LT(angularError(filter.orientation(), Quaternion::Identity()), 1e-15);
}
