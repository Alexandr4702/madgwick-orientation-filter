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

- C++17 compiler
- CMake 3.16 or newer
- Eigen 3
- GoogleTest when tests are enabled

CMake downloads Eigen when it is not available locally. GoogleTest is fetched
only when `MADGWICK_BUILD_TESTS` is enabled.

## Coordinate convention

The filter stores an internal unit quaternion

$$
q = q_{BW},
$$

and returns its inverse:

$$
Q = q^{-1}.
$$

`get_orientation()` returns `Q`.

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

Reference and measured vectors must be normalized before calling `update()`.
Angular velocity must be expressed in radians per second, and `dt` in seconds.

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
\dot q_\omega = \frac{1}{2}q \otimes \Omega,
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

For

$$
q=[w,x,y,z]^T, \qquad d=[d_x,d_y,d_z]^T,
$$

and quaternion columns ordered as $[w,x,y,z]$, the `Jacobian()` helper in the
current source evaluates

$$
J(q,d)=2
\begin{bmatrix}
d_yz-d_zy & d_yy+d_zz & -2d_xy+d_yx-d_zw & -2d_xz+d_yw+d_zx \\
-d_xz+d_zx & d_xy-2d_yx+d_zw & d_xx+d_zz & -d_xw-2d_yz+d_zy \\
d_xy-d_yx & d_xz-d_yw-2d_zx & d_xw+d_yz-2d_zy & d_xx+d_yy
\end{bmatrix}.
$$

Written element by element, this is

$$
\begin{aligned}
J_{11}&=2(d_yz-d_zy), &
J_{12}&=2(d_yy+d_zz), \\
J_{13}&=2(-2d_xy+d_yx-d_zw), &
J_{14}&=2(-2d_xz+d_yw+d_zx), \\
J_{21}&=2(-d_xz+d_zx), &
J_{22}&=2(d_xy-2d_yx+d_zw), \\
J_{23}&=2(d_xx+d_zz), &
J_{24}&=2(-d_xw-2d_yz+d_zy), \\
J_{31}&=2(d_xy-d_yx), &
J_{32}&=2(d_xz-d_yw-2d_zx), \\
J_{33}&=2(d_xw+d_yz-2d_zy), &
J_{34}&=2(d_xx+d_yy).
\end{aligned}
$$

Products such as $d_yz$ mean ordinary scalar multiplication. The measured
vector $s$ does not appear in $J$ because it is constant with respect to the
quaternion; it enters through the residual $f$.

The legacy `update()` path does not call `Jacobian()` directly. It calls
`grad_function()`, which contains a symbolically expanded gradient. The matrix
above documents the existing helper exactly; the simulation and
characterization tests protect the behavior of the expanded update path.

The source contains the expanded analytic form of this gradient, avoiding a
numerical Jacobian at runtime.

For two observations, the gradients are added:

$$
g = g_1 + g_2.
$$

The current implementation intentionally preserves its original behavior and
does not normalize this combined gradient. Consequently, the effective
correction strength depends on the residual magnitude and the geometry of the
reference vectors.

### Corrected state equation

The complete derivative is

$$
\dot q = \frac{1}{2}q\otimes\Omega - \beta g,
$$

where $\beta$ controls the strength of the vector correction. The current
default is

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
\operatorname{vec}\left(2q^{-1}\otimes g\right).
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

Madgwick_filter filter;

Vec3 gravity_world(0.0, 0.0, 1.0);
Vec3 magnetic_world(1.0, 0.0, 0.0);

Vec3 gravity_body = accelerometer.normalized();
Vec3 magnetic_body = magnetometer.normalized();
Vec3 gyro_rad_s = gyroscope;

filter.update(
    gravity_world,
    magnetic_world,
    gravity_body,
    magnetic_body,
    gyro_rad_s,
    dt_seconds);

Quat orientation = filter.get_orientation();
Vec3 estimated_bias = filter.getOmega_bias();
```

### One-vector update

Use this mode when only one absolute direction is available:

```cpp
filter.update(
    gravity_world,
    accelerometer.normalized(),
    gyroscope_rad_s,
    dt_seconds);
```

The unconstrained rotation about the supplied direction is propagated only by
the gyroscope and will therefore drift over time.

### Initial gyroscope bias

If a calibration value is available, it can be supplied before filtering:

```cpp
filter.setOmega_bias(calibrated_bias_rad_s);
```

### Angular velocity from quaternion samples

```cpp
Vec3 omega = Madgwick_filter::get_omega_from_quat(
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

To use only system-installed dependencies:

```sh
cmake -S . -B build \
  -DMADGWICK_BUILD_TESTS=ON \
  -DMADGWICK_FETCH_DEPENDENCIES=OFF
```

Projects using CMake can include the library with:

```cmake
add_subdirectory(path/to/madgwick-orientation-filter)
target_link_libraries(your_target PRIVATE Madgwick::madgwick)
```

## Tests

The GoogleTest suite contains both characterization tests and independent
sensor simulations.

### Characterization tests

These tests lock down the current numerical behavior so that future refactors
can be checked for unintended changes:

- identity preservation with zero angular velocity;
- the exact existing Euler-integration result for one gyroscope step;
- bit-for-bit repeatability of a fixed reference scenario;
- angular velocity recovered from a quaternion delta;
- exact round-trip behavior of the gyroscope-bias accessors.

### Simulation tests

The simulations integrate an independent ground-truth quaternion, synthesize
ideal vector observations from it, feed those observations to the filter, and
compare the estimated attitude with the known truth.

The suite covers:

- constant three-axis rotation with ideal gyroscope and direction sensors;
- convergence from an initial attitude error using two fixed directions;
- stationary operation with a constant gyroscope bias.

The ground-truth propagation uses `Eigen::AngleAxisd`, while the filter uses
its original first-order quaternion integration. This separation prevents the
simulation from merely duplicating the implementation under test.

## Current limitations

- Input vectors are expected to be normalized by the caller.
- The implementation does not currently reject zero vectors, non-finite
  values, or non-positive time steps.
- The correction gain and bias gain are fixed in the class implementation.
- Euler integration is accurate only when the sampling interval is small
  relative to angular velocity.
- Accelerometer measurements are reliable gravity observations only when
  linear acceleration is sufficiently small.
- Magnetometer measurements require calibration and can be disturbed by nearby
  magnetic materials and currents.
