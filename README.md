# Generalized Madgwick Orientation Filter

A compact C++ implementation of a Madgwick-style attitude estimator for
arbitrary reference vectors.

Traditional Madgwick implementations are usually hard-coded for gravity and
the Earth's magnetic field. This implementation accepts any known direction,
including:

- gravity measured by an accelerometer;
- the Earth's magnetic field;
- a Sun direction from a Sun sensor;
- a star-tracker line of sight;
- a visual landmark direction;
- any other direction known in both the reference and body frames.

The filter supports either one reference vector or two reference vectors and
includes online gyroscope-bias estimation.

## Dependencies

- C++20 compiler
- CMake 3.16 or newer
- Eigen 3 Git submodule
- System-installed GoogleTest when tests are enabled

Eigen is included as the `third_party/eigen` Git submodule.

## Coordinate convention

The filter stores an internal unit quaternion

$$
q = q_{BW},
$$

and returns its inverse:

$$
Q = q^{-1}.
$$

`orientation()` returns `Q`.

Here, $q_{BW}$ rotates a world-frame vector into the body frame, while the
returned $Q=q_{WB}$ rotates a body-frame vector into the world frame. For a
reference direction $d$ and its sensor measurement $s$, the convention is

$$
s = Q^{-1}\,d\,Q = q\,d\,q^{-1}.
$$

Vectors are represented as pure quaternions when they participate in
quaternion products:

$$
d_q = [0, d_x, d_y, d_z]^T.
$$

Reference and measured vectors are normalized by `update()`. Zero and non-finite
vectors are rejected. Angular velocity must be in radians per second and `dt`
must be finite and positive.

## Filter model

### Gyroscope propagation

Given the measured body angular velocity $\omega_m$ and estimated bias $b$,
the corrected angular velocity is

$$
\omega = \omega_m - b.
$$

Using the pure quaternion

$$
\Omega = [0, \omega_x, \omega_y, \omega_z]^T,
$$

the gyroscope contribution to the quaternion derivative is

$$
\dot q_\omega = -\frac{1}{2}\Omega \otimes q,
$$

where $\otimes$ denotes quaternion multiplication.

### Vector observation error

For a known reference vector $d$ and measured vector $s$, the observation
residual can be written as

$$
f(q,d,s) = \operatorname{vec}\left(q\otimes d_q\otimes q^{-1}\right)-s.
$$

The correction direction is obtained from the gradient of the squared
observation error:

$$
g(q) = J(q,d)^T f(q,d,s),
$$

with

$$
J(q,d)=\frac{\partial f(q,d,s)}{\partial q}.
$$

The implementation evaluates this analytic Jacobian directly. A randomized
test compares all four gradient components against central finite differences.

For two observations, the gradients are added:

$$
g = g_1 + g_2.
$$

The current implementation intentionally preserves its original behavior and
does not normalize this combined gradient. Consequently, the effective
correction strength depends on the residual magnitude and the geometry of the
reference vectors.

### Corrected state equation

Expressed in this quaternion convention, the original Madgwick
normalized-gradient form is:

$$
\dot q = -\frac{1}{2}\Omega\otimes q
- \beta\frac{g}{\lVert g\rVert}.
$$

This generalized implementation uses the raw gradient instead:

$$
\dot q = -\frac{1}{2}\Omega\otimes q - \beta g,
$$

so the correction strength depends on the measurement error and reference
vector geometry. This gives a smooth, vanishing correction near the solution
and avoids division by a near-zero norm, but makes convergence and tuning
depend on the number, scaling, and geometry of observations. The normalized
form has a correction magnitude fixed by $\beta$, but may amplify sensor noise
near the solution. The current default is

$$
\beta = 0.1.
$$

The implementation uses a first-order Euler step followed by normalization:

$$
q_{k+1} = \operatorname{normalize}\left(q_k + \dot q_k\Delta t\right).
$$

Normalization keeps the quaternion on the unit sphere despite integration and
floating-point error.

### Gyroscope-bias estimation

The angular correction associated with the gradient is calculated as

$$
\omega_\varepsilon =
-\operatorname{vec}\left(2g\otimes q^{-1}\right).
$$

The bias estimate is then updated using

$$
b_{k+1} = b_k + k_b\,\omega_\varepsilon\Delta t,
$$

where the implementation currently uses

$$
k_b = 0.001.
$$

## Observability

A single direction provides only two independent attitude constraints.
Rotation about that direction is unobservable. For example, gravity alone can
estimate roll and pitch, but not heading.

Two non-collinear directions provide full attitude observability. Nearly
parallel vectors should be avoided because they produce weak numerical
conditioning and poor heading accuracy.

## Usage

### Two-vector update

This is the preferred mode when full attitude estimation is required.

```cpp
#include "Madgwick.h"

madgwick::Filter filter;

madgwick::Vector3 gyro_rad_s = gyroscope;

filter.update(
    {madgwick::Vector3(0.0, 0.0, 1.0), accelerometer},
    {madgwick::Vector3(1.0, 0.0, 0.0), magnetometer},
    gyro_rad_s,
    dt_seconds);

madgwick::Quaternion orientation = filter.orientation();
madgwick::Vector3 estimated_bias = filter.gyro_bias();
```

### One-vector update

Use this mode when only one absolute direction is available:

```cpp
filter.update(
    {gravity_world, accelerometer},
    gyroscope_rad_s,
    dt_seconds);
```

The unconstrained rotation about the supplied direction is propagated only by
the gyroscope and will therefore drift over time.

### Initial gyroscope bias

If a calibration value is available, it can be supplied before filtering:

```cpp
filter.set_gyro_bias(calibrated_bias_rad_s);
```

### Angular velocity from quaternion samples

```cpp
madgwick::Vector3 omega = madgwick::Filter::angular_velocity(
    current_orientation,
    previous_orientation,
    dt_seconds);
```

## Use cases

### IMU attitude estimation

Use gravity and magnetic field directions with accelerometer, magnetometer,
and gyroscope measurements to estimate a complete 3D orientation.

### Spacecraft attitude determination

Use a Sun vector together with a magnetic-field vector or star-tracker vector.
The generalized observation interface avoids assumptions specific to
terrestrial IMUs.

### Robotics and visual navigation

Fuse angular velocity with known landmark, surface-normal, gravity, or heading
directions. Measurements may come from different sensors as long as they are
expressed in a consistent body frame.

### Heading-free stabilization

Use the single-vector overload with gravity when only roll and pitch
stabilization are required and yaw drift is acceptable.

## Build

```sh
cmake -S . -B build -DMADGWICK_BUILD_TESTS=ON
cmake --build build --config Release
```

Run the test suite with:

```sh
ctest --test-dir build -C Release --output-on-failure
```

Clone with the Eigen submodule:

```sh
git clone --recurse-submodules <repository-url>
```

For an existing checkout:

```sh
git submodule update --init --recursive
```

Projects using CMake can include the library with:

```cmake
add_subdirectory(path/to/madgwick-orientation-filter)
target_link_libraries(your_target PRIVATE madgwick)
```

## Tests

The GoogleTest suite checks numerical behavior and runs deterministic sensor
simulations with motion, noise, timing jitter, and gyroscope bias.

## Current limitations

- Euler integration is accurate only when the sampling interval is small
  relative to angular velocity.
- Accelerometer measurements are reliable gravity observations only when
  linear acceleration is sufficiently small.
- Magnetometer measurements require calibration and can be disturbed by nearby
  magnetic materials and currents.
