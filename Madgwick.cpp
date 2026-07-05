//============================================================================
// Name        : Madgwick_test.cpp
// Author      : 
// Version     :
// Copyright   : Your copyright notice
// Description : Hello World in C++, Ansi-style
//============================================================================

#include "Madgwick.h"

/*
 *@brief: function for work with quaternions
 *
*/
//----------------------------------------------------------------------------------------------------------------
Quat operator* (Quat& q1, double Scalar )
{
    return Quat(q1.w() * Scalar, q1.x() * Scalar, q1.y() * Scalar, q1.z() * Scalar);
}
Quat operator* (Quat&& q1, double Scalar )
{
    return Quat(q1.w() * Scalar, q1.x() * Scalar, q1.y() * Scalar, q1.z() * Scalar);
}
Quat operator* (double Scalar, Quat& q1 ) {
    return Quat(q1.w() * Scalar, q1.x() * Scalar, q1.y() * Scalar, q1.z() * Scalar);
}
Quat operator* (double Scalar, Quat&& q1 ) {
    return Quat(q1.w() * Scalar, q1.x() * Scalar, q1.y() * Scalar, q1.z() * Scalar);
}
//----------------------------------------------------------------------------------------------------------------
Quat operator/ (Quat& q1, double Scalar )
{
    return Quat(q1.w() / Scalar, q1.x() / Scalar, q1.y() / Scalar, q1.z() / Scalar);
}
Quat operator/ (Quat&& q1, double Scalar )
{
    return Quat(q1.w() / Scalar, q1.x() / Scalar, q1.y() / Scalar, q1.z() / Scalar);
}
Quat operator/ (double Scalar, Quat& q1 ) {
    return Quat(q1.w() / Scalar, q1.x() / Scalar, q1.y() / Scalar, q1.z() / Scalar);
}
Quat operator/ (double Scalar, Quat&& q1 ) {
    return Quat(q1.w() / Scalar, q1.x() / Scalar, q1.y() / Scalar, q1.z() / Scalar);
}
//----------------------------------------------------------------------------------------------------------------???
Quat operator+ (Quat& q1, double Scalar )
{
    return Quat(q1.w() + Scalar, q1.x() + Scalar, q1.y() + Scalar, q1.z() + Scalar);
}
Quat operator+ (Quat&& q1, double Scalar )
{
    return Quat(q1.w() + Scalar, q1.x() + Scalar, q1.y() + Scalar, q1.z() + Scalar);
}
Quat operator+ (double Scalar, Quat& q1 )
{
    return Quat(q1.w() + Scalar, q1.x() + Scalar, q1.y() + Scalar, q1.z() + Scalar);
}
Quat operator+ (double Scalar, Quat&& q1 )
{
    return Quat(q1.w() + Scalar, q1.x() + Scalar, q1.y() + Scalar, q1.z() + Scalar);
}

Quat operator+ (Quat& q1, Quat& q2 )
{
    return Quat(q1.w() + q2.w(), q1.x() + q2.x(), q1.y() + q2.y(), q1.z() + q2.z());
}
Quat operator+ (Quat&& q1, Quat& q2 )
{
    return Quat(q1.w() + q2.w(), q1.x() + q2.x(), q1.y() + q2.y(), q1.z() + q2.z());
}
Quat operator+ (Quat& q1, Quat&& q2 )
{
    return Quat(q1.w() + q2.w(), q1.x() + q2.x(), q1.y() + q2.y(), q1.z() + q2.z());
}
Quat operator+ (Quat&& q1, Quat&& q2 )
{
    return Quat(q1.w() + q2.w(), q1.x() + q2.x(), q1.y() + q2.y(), q1.z() + q2.z());
}
//----------------------------------------------------------------------------------------------------------------
Quat operator- (Quat& q1, double Scalar )
{
    return Quat(q1.w() - Scalar, q1.x() - Scalar, q1.y() - Scalar, q1.z() - Scalar);
}
Quat operator- (Quat&& q1, double Scalar )
{
    return Quat(q1.w() - Scalar, q1.x() - Scalar, q1.y() - Scalar, q1.z() - Scalar);
}
Quat operator- (double Scalar, Quat& q1 )
{
    return Quat(q1.w() - Scalar, q1.x() - Scalar, q1.y() - Scalar, q1.z() - Scalar);
}
Quat operator- (double Scalar, Quat&& q1 )
{
    return Quat(q1.w() - Scalar, q1.x() - Scalar, q1.y() - Scalar, q1.z() - Scalar);
}

Quat operator- (Quat& q1, Quat& q2 )
{
    return Quat(q1.w() - q2.w(), q1.x() - q2.x(), q1.y() - q2.y(), q1.z() - q2.z());
}
Quat operator- (Quat&& q1, Quat& q2 )
{
    return Quat(q1.w() - q2.w(), q1.x() - q2.x(), q1.y() - q2.y(), q1.z() - q2.z());
}
Quat operator- (Quat& q1, Quat&& q2 )
{
    return Quat(q1.w() - q2.w(), q1.x() - q2.x(), q1.y() - q2.y(), q1.z() - q2.z());
}
Quat operator- (Quat&& q1, Quat&& q2 )
{
    return Quat(q1.w() - q2.w(), q1.x() - q2.x(), q1.y() - q2.y(), q1.z() - q2.z());
}

void Madgwick_filter::update(Vec3 sun_ref, Vec3 magn_ref, Vec3 sun_meas, Vec3 magn_meas, Vec3 omega, double dt)
{
    Quat grad_sun = grad_function(sun_ref, sun_meas, q_1);
    Quat grad_magn = grad_function(magn_ref, magn_meas, q_1);
    Quat grad = grad_sun + grad_magn;
    //       grad.normalize();
    Vec3 omega_eps_t = (2 * q_1.inverse() * grad).vec();
    omega_bias += omega_eps_t * dt * 0.001;
    Quat q_dervative_avs = q_dervative_omega(q_1, omega - omega_bias );

    Quat q_dervative;
    q_dervative = q_dervative_avs - betta * grad ;
    q_1 =  q_1 + q_dervative * dt;
    q_1.normalize();
}

void Madgwick_filter::update(Vec3 sun_ref, Vec3 sun_meas, Vec3 omega, double dt)
{
    Quat grad_sun = grad_function(sun_ref, sun_meas, q_1);
    Quat grad = grad_sun;
    //       grad.normalize();
    Vec3 omega_eps_t = (2 * q_1.inverse() * grad).vec();
    omega_bias += omega_eps_t * dt * 0.001;
    Quat q_dervative_avs = q_dervative_omega(q_1, omega - omega_bias );

    Quat q_dervative;
    q_dervative = q_dervative_avs - betta * grad ;
    q_1 =  q_1 + q_dervative * dt;
    q_1.normalize();
}

Quat Madgwick_filter::get_orientation()
{
    return q_1.inverse();
}

Vec3 Madgwick_filter::getOmega_bias() const
{
return omega_bias;
}

void Madgwick_filter::setOmega_bias(const Vec3 &value)
{
omega_bias = value;
}

Vec3 Madgwick_filter::get_omega_from_quat(Quat q, Quat q_1, double dt)
{
    Quat delata_q = q_1.inverse() * q;
    return acos(delata_q.w()) * 2 * delata_q.vec().normalized() / dt;
    return (2 * (q - q_1) * q.conjugate() / dt).vec();
}

Quat Madgwick_filter::omega_integration (Quat q_1,Vec3 omega, double dt)
{
    Quat q_dervative = q_dervative_omega(q_1, omega);
    Quat ret = q_1 + q_dervative * dt;
    return ret;
}

Quat Madgwick_filter::func(Quat q, Vec3 d, Vec3 s)
 {
//        Quat q;
     Vec3 ref = d ;
     Vec3 frame = s;
     Quat f;
//        double q1q2 = q.w() * q.x();
//        double q1q3 = q.w() * q.y();
//        double q1q4 = q.w() * q.z();

//        double q2q2 = q.x() * q.x();
//        double q2q3 = q.x() * q.y();
//        double q2q4 = q.x() * q.z();

//        double q3q3 = q.y() * q.y();
//        double q3q4 = q.y() * q.z();

//        double q4q4 = q.z() * q.z();

//        double _2dx = 2.0 * ref.x();
//        double _2dy = 2.0 * ref.y();
//        double _2dz = 2.0 * ref.z();

//        f.w() = 0;

//        f.x() =
//                _2dx * (0.5 - q3q3 - q4q4) +
//                _2dy * (q1q4 + q2q3) +
//                _2dz * (q2q4 - q1q3)
//                - frame.x();
//        f.y() =
//                _2dx * (q2q3 - q1q4) +
//                _2dy * (0.5 - q2q2 - q4q4) +
//                _2dz * (q1q2 + q3q4)
//                - frame.y();
//        f.z() =
//                _2dx * (q1q3 + q2q4) +
//                _2dy * (q3q4 - q1q2) +
//                _2dz * (0.5 - q2q2 - q3q3)
//                - frame.z();
//        return f;

     Quat ref_q(0,ref.x (),ref.y (),ref.z ());
     Quat frame_q(0,frame.x (),frame.y (),frame.z ());
     return q.conjugate () * ref_q * q - frame_q;
 }

Eigen::Matrix<double, 3, 4> Madgwick_filter::Jacobian(Quat q, Vec3 d)
{
//        Vec3 d;
//        Quat q;
    double q1 = q.w();
    double q2 = q.x();
    double q3 = q.y();
    double q4 = q.z();

    Eigen::Matrix<double, 3, 4> ret;
    ret << ( d.y() * q4 - d.z() * q3), (d.y() * q3 + d.z() * q4)                 , ( - 2.0 * d.x() * q3 + d.y() * q2 - d.z() * q1), (-2 * d.x() * q4 + d.y() * q1 + d.z() * q2),
           (-d.x() * q4 + d.z() * q2), (d.x() * q3 - 2 * d.y() * q2 + d.z() * q1), d.x() * q2 + d.z() * q4                        , (-d.x() * q1 - 2  * d.y() * q4 + d.z() * q3),
             d.x() * q3 - d.y() * q2 ,  d.x() * q4 - d.y() * q1 - 2 * d.z() * q2 , d.x() * q1 + d.y() * q4 - 2 * d.z() * q3       ,  d.x() * q2 + d.y() * q3;
    ret = ret * 2;
    return ret;
}

Quat Madgwick_filter::grad_function(Vec3 ref, Vec3 meas, Quat q)
{
    Quat ret;
    double q_module = (q * q.conjugate()).w();
    double measx_2 = meas.x() * meas.x();
    double measy_2 = meas.y() * meas.y();
    ret.w() = 2 * ( q.w() * (measx_2 * q_module + measy_2 * q_module) + meas.z() * (meas.z() * q.w() * q_module + q.y()*ref.x() - q.x() * ref.y() - q.w() * ref.z()) - meas.y() * (q.z() * ref.x() + q.w() * ref.y() - q.x() * ref.z()) - meas.x() * (q.w() * ref.x() - q.z() * ref.y() + q.y() * ref.z()) );
    ret.x() = 2 * (measx_2 * q.x() * q_module +
               measy_2 * q.x() * q_module +
               meas.y() * (-q.y() * ref.x() + q.x() * ref.y() + q.w() * ref.z()) +
               meas.z() * (meas.z() * q.x() * q_module - q.z() * ref.x() - q.w() * ref.y() +
                  q.x() * ref.z()) - meas.x() * (q.x() * ref.x() + q.y() * ref.y() + q.z() * ref.z()));
    ret.y() = 2 * (measx_2 * q.y() * q_module +
               measy_2 * q.y() * q_module +
               meas.x() * (q.y() * ref.x() - q.x() * ref.y() - q.w() * ref.z()) +
               meas.z() * (meas.z() * q.y() * q_module + q.w() * ref.x() - q.z() * ref.y() +
                  q.y() * ref.z()) - meas.y() * (q.x() * ref.x() + q.y() * ref.y() + q.z() * ref.z()));
    ret.z() = 2 * (measx_2 * q.z() * q_module +
               measy_2 * q.z() * q_module +
               meas.x() * (q.z() * ref.x() + q.w() * ref.y() - q.x() * ref.z()) -
               meas.y() * (q.w() * ref.x() - q.z() * ref.y() + q.y() * ref.z()) +
               meas.z() * (meas.z() * q.z() * q_module - q.x() * ref.x() - q.y() * ref.y() -
                  q.z() * ref.z()));
    return ret;
}

Quat Madgwick_filter::q_dervative_omega (Quat q_1, Vec3 omega)
{
    Quat omega_q;
    Quat ret;
    omega_q.coeffs().segment(0,3) = omega;
    omega_q.w() = 0;
    Quat q_dervative = 0.5 * q_1 *omega_q;
    return q_dervative;
}
