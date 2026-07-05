#include "Madgwick.h"

#include <gtest/gtest.h>
#include <algorithm>
#include <array>
#include <cmath>
#include <random>

namespace {
double angularError(const Quat& actual, const Quat& expected)
{
    const double dot = std::min(1.0, std::abs(actual.normalized().dot(expected.normalized())));
    return 2.0 * std::acos(dot);
}

Quat integrateTruth(const Quat& orientation, const Vec3& omega, double dt)
{
    const double angle = omega.norm() * dt;
    if (angle == 0.0)
        return orientation;
    return (Quat(Eigen::AngleAxisd(angle, omega.normalized())) * orientation).normalized();
}

Vec3 noisyDirection(const Vec3& direction, std::mt19937& generator,
                    std::normal_distribution<double>& noise)
{
    return (direction + Vec3(noise(generator), noise(generator), noise(generator))).normalized();
}

std::array<double, 7> runReferenceScenario()
{
    Madgwick_filter filter;
    const Vec3 ref1(1.0, 0.0, 0.0);
    const Vec3 ref2 = Vec3(0.5, 0.0, 0.5).normalized();
    Quat orientation(1.0, 0.0, 0.0, 0.0);
    const Vec3 omega(0.01, -0.02, 0.03);

    for (int i = 0; i < 100; ++i) {
        const Vec3 meas1 = orientation * ref1;
        const Vec3 meas2 = orientation * ref2;
        filter.update(ref1, ref2, meas1, meas2, omega, 0.01);

        Quat omega_step(1.0, 0.005 * omega.x(), 0.005 * omega.y(),
                        0.005 * omega.z());
        orientation = (orientation * omega_step).normalized();
    }

    const Quat q = filter.get_orientation();
    const Vec3 bias = filter.getOmega_bias();
    return {q.w(), q.x(), q.y(), q.z(), bias.x(), bias.y(), bias.z()};
}
} // namespace

TEST(Characterization, IdentityIsPreservedWithZeroAngularVelocity)
{
    Madgwick_filter filter;
    for (int i = 0; i < 100; ++i)
        filter.update(Vec3::UnitX(), Vec3::UnitY(), Vec3::UnitX(),
                      Vec3::UnitY(), Vec3::Zero(), 0.01);

    const Quat q = filter.get_orientation();
    EXPECT_DOUBLE_EQ(q.w(), 1.0);
    EXPECT_DOUBLE_EQ(q.x(), 0.0);
    EXPECT_DOUBLE_EQ(q.y(), 0.0);
    EXPECT_DOUBLE_EQ(q.z(), 0.0);
}

TEST(Characterization, OneGyroStepMatchesExistingEulerIntegration)
{
    Madgwick_filter filter;
    filter.update(Vec3::UnitZ(), Vec3::UnitZ(), Vec3(1.0, 2.0, 3.0), 0.1);

    // Existing implementation integrates (1, 0.05, 0.10, 0.15), normalizes
    // it, and returns its inverse. These equations intentionally mirror that
    // behavior rather than prescribing a new integration method.
    const double norm = std::sqrt(1.0 + 0.05*0.05 + 0.10*0.10 + 0.15*0.15);
    const Quat q = filter.get_orientation();
    EXPECT_DOUBLE_EQ(q.w(), 1.0 / norm);
    EXPECT_DOUBLE_EQ(q.x(), -0.05 / norm);
    EXPECT_DOUBLE_EQ(q.y(), -0.10 / norm);
    EXPECT_DOUBLE_EQ(q.z(), -0.15 / norm);
}

TEST(Characterization, ReferenceScenarioIsBitwiseRepeatable)
{
    const auto first = runReferenceScenario();
    const auto second = runReferenceScenario();
    for (std::size_t i = 0; i < first.size(); ++i)
        EXPECT_DOUBLE_EQ(first[i], second[i]) << "component " << i;
}

TEST(Simulation, TracksConstantRotationWithIdealSensors)
{
    Madgwick_filter filter;
    Quat truth = Quat::Identity();
    const Vec3 omega(0.08, -0.12, 0.20);
    const Vec3 ref1 = Vec3::UnitZ();
    const Vec3 ref2 = Vec3(0.8, 0.2, 0.1).normalized();
    constexpr double dt = 0.002;

    for (int i = 0; i < 5000; ++i) {
        // get_orientation() returns the inverse of the internally integrated
        // quaternion, hence its positive gyro convention is -omega.
        truth = integrateTruth(truth, -omega, dt);
        const Vec3 meas1 = truth.conjugate() * ref1;
        const Vec3 meas2 = truth.conjugate() * ref2;
        filter.update(ref1, ref2, meas1, meas2, omega, dt);
    }

    EXPECT_LT(angularError(filter.get_orientation(), truth), 5e-4);
}

TEST(Simulation, CorrectsAnInitialAttitudeError)
{
    Madgwick_filter filter;
    const Quat truth(Eigen::AngleAxisd(0.35, Vec3(1.0, -2.0, 0.5).normalized()));
    const Vec3 ref1 = Vec3::UnitZ();
    const Vec3 ref2 = Vec3::UnitX();
    const Vec3 meas1 = truth.conjugate() * ref1;
    const Vec3 meas2 = truth.conjugate() * ref2;

    for (int i = 0; i < 5000; ++i)
        filter.update(ref1, ref2, meas1, meas2, Vec3::Zero(), 0.002);

    // The legacy, unnormalized gradient settles to a small residual error.
    EXPECT_LT(angularError(filter.get_orientation(), truth), 5e-3);
}

TEST(Simulation, EstimatesConstantGyroscopeBiasWhileStationary)
{
    Madgwick_filter filter;
    const Vec3 gyro_bias(0.004, -0.003, 0.002);
    const Vec3 ref1 = Vec3::UnitZ();
    const Vec3 ref2 = Vec3::UnitX();

    for (int i = 0; i < 30000; ++i)
        filter.update(ref1, ref2, ref1, ref2, gyro_bias, 0.002);

    EXPECT_LT(angularError(filter.get_orientation(), Quat::Identity()), 1e-2);
    EXPECT_LT((filter.getOmega_bias() - gyro_bias).norm(), gyro_bias.norm());
}

TEST(Simulation, TracksTimeVaryingThreeAxisMotionWithSamplePeriodJitter)
{
    Madgwick_filter filter;
    Quat truth = Quat::Identity();
    const Vec3 ref1 = Vec3(0.2, -0.3, 0.93).normalized();
    const Vec3 ref2 = Vec3(-0.7, 0.6, 0.2).normalized();
    double time = 0.0;

    for (int i = 0; i < 12000; ++i) {
        const double dt = 0.0015 + 0.0005 * (1.0 + std::sin(0.017 * i));
        const double middle = time + 0.5 * dt;
        const Vec3 omega(0.45 * std::sin(0.71 * middle),
                         -0.35 * std::cos(0.43 * middle),
                         0.25 + 0.20 * std::sin(0.19 * middle));

        truth = integrateTruth(truth, -omega, dt);
        filter.update(ref1, ref2, truth.conjugate() * ref1,
                      truth.conjugate() * ref2, omega, dt);
        time += dt;
    }

    EXPECT_LT(angularError(filter.get_orientation(), truth), 2e-3);
}

TEST(Simulation, RemainsAccurateWithGyroscopeAndDirectionNoise)
{
    Madgwick_filter filter;
    Quat truth = Quat::Identity();
    const Vec3 ref1 = Vec3::UnitZ();
    const Vec3 ref2 = Vec3(0.8, -0.3, 0.15).normalized();
    std::mt19937 generator(0x4d414447);
    std::normal_distribution<double> direction_noise(0.0, 0.004);
    std::normal_distribution<double> gyro_noise(0.0, 0.002);
    constexpr double dt = 0.002;

    for (int i = 0; i < 15000; ++i) {
        const double time = i * dt;
        const Vec3 omega(0.18 * std::sin(0.37 * time),
                         0.12 * std::cos(0.23 * time),
                         -0.16 * std::sin(0.51 * time));
        truth = integrateTruth(truth, -omega, dt);

        const Vec3 gyro = omega + Vec3(gyro_noise(generator),
                                       gyro_noise(generator),
                                       gyro_noise(generator));
        const Vec3 meas1 = noisyDirection(truth.conjugate() * ref1, generator,
                                          direction_noise);
        const Vec3 meas2 = noisyDirection(truth.conjugate() * ref2, generator,
                                          direction_noise);
        filter.update(ref1, ref2, meas1, meas2, gyro, dt);
    }

    EXPECT_LT(angularError(filter.get_orientation(), truth), 1e-2);
}

TEST(Simulation, TracksManeuversWhileEstimatingConstantGyroscopeBias)
{
    Madgwick_filter filter;
    Quat truth = Quat::Identity();
    const Vec3 ref1 = Vec3(0.3, 0.1, 0.95).normalized();
    const Vec3 ref2 = Vec3(-0.2, 0.9, 0.35).normalized();
    const Vec3 gyro_bias(0.006, -0.004, 0.003);
    constexpr double dt = 0.002;

    for (int i = 0; i < 30000; ++i) {
        const double time = i * dt;
        const Vec3 omega(0.20 * std::sin(0.29 * time),
                         -0.15 * std::sin(0.47 * time),
                         0.10 * std::cos(0.31 * time));
        truth = integrateTruth(truth, -omega, dt);
        filter.update(ref1, ref2, truth.conjugate() * ref1,
                      truth.conjugate() * ref2, omega + gyro_bias, dt);
    }

    EXPECT_LT(angularError(filter.get_orientation(), truth), 2e-2);
    EXPECT_LT((filter.getOmega_bias() - gyro_bias).norm(), gyro_bias.norm());
}

TEST(Simulation, ConvergesFromLargeErrorsAroundDifferentAxes)
{
    const Vec3 ref1 = Vec3(0.1, -0.2, 0.97).normalized();
    const Vec3 ref2 = Vec3(0.8, 0.5, -0.1).normalized();
    const std::array<Vec3, 4> axes = {
        Vec3::UnitX(), Vec3::UnitY(), Vec3::UnitZ(), Vec3(1.0, -2.0, 3.0).normalized()
    };

    for (const Vec3& axis : axes) {
        Madgwick_filter filter;
        const Quat truth(Eigen::AngleAxisd(1.4, axis));
        const Vec3 meas1 = truth.conjugate() * ref1;
        const Vec3 meas2 = truth.conjugate() * ref2;

        for (int i = 0; i < 20000; ++i)
            filter.update(ref1, ref2, meas1, meas2, Vec3::Zero(), 0.002);

        EXPECT_LT(angularError(filter.get_orientation(), truth), 3e-2)
            << "axis: " << axis.transpose();
    }
}

TEST(Characterization, QuaternionDeltaProducesExpectedAngularVelocity)
{
    const double angle = 0.2;
    const Quat previous = Quat::Identity();
    const Quat current(Eigen::AngleAxisd(angle, Vec3::UnitY()));
    const Vec3 omega = Madgwick_filter::get_omega_from_quat(current, previous, 0.1);
    EXPECT_DOUBLE_EQ(omega.x(), 0.0);
    EXPECT_DOUBLE_EQ(omega.y(), 1.9999999999999889);
    EXPECT_DOUBLE_EQ(omega.z(), 0.0);
}

TEST(Characterization, LegacyBiasAccessorsRoundTripExactly)
{
    Madgwick_filter filter;
    const Vec3 bias(0.125, -0.25, 0.5);
    filter.setOmega_bias(bias);
    EXPECT_TRUE(filter.getOmega_bias() == bias);
}
